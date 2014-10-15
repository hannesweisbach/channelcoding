#include <algorithm>
#include <sstream>
#include <numeric>
#include <iomanip>
#include <cstdlib>
#include <clocale>
#include <typeinfo>
#include <stdexcept>

#include "rs.h"
#include "cyclic.h"
#include "center.h"

rs::rs(unsigned q, uint64_t modular_polynomial, unsigned dmin)
    : field(q, modular_polynomial), n((1 << q) - 1), k(dmin - 1),
      l(n - dmin + 1), dmin(dmin), mu(1), fk((dmin - 1) / 2),
      generator(field, { 1 }) {
  for (int i = mu; i < mu + dmin - 1; i++) {
    gf_polynomial root(field, { i, 0 }, true);
    generator *= root;
  }
  std::cout << "(" << n << ", " << l << ", " << dmin << ")- RS code over GF(2^"
            << q << ") with generator " << generator << std::endl;
}

gf_polynomial rs::pzg(const std::vector<gf_element> &syndromes) const {
  std::cout << std::endl << "Calculating Σ(x) with PZG:" << std::endl;
  for (unsigned v = fk; v; v--) {
    gf_matrix eq_system(field);
    for (auto it = std::cbegin(syndromes); it != std::cbegin(syndromes) + v;
         ++it) {
      gf_polynomial poly(field);
      for (auto cp_it = it; cp_it != it + v + 1; ++cp_it)
        poly.push_front(*cp_it);
      eq_system.push_back(poly);
    }

    std::cout << "Linear equation system to be solved: " << std::endl;
    std::cout << eq_system << std::endl;

    // TODO Check if the determinant is non-zero.
    try {
      std::cout << "Reduced row echelon form of the linear equation system:"
                << std::endl;

      std::cout << eq_system.reduced_echelon_form() << std::endl;

      std::cout << "Solution vector of the linear equation system:"
                << std::endl;
      auto solution = eq_system.solution();
      std::cout << solution << std::endl << std::endl;

      return solution;
    }
    catch (const std::exception &e) {
      std::cout << "Error for v = " << v << ": " << typeid(e).name() << " "
                << e.what() << std::endl;
    }
  }

  /* if v = 1, check all equations for the same solution for sigma_1 */
  auto sigma = syndromes.front();
  const auto one = field.power_to_polynomial(0);
  for (auto it = std::begin(syndromes); it != std::end(syndromes) - 1; ++it) {
    auto sigma_ = *(it + 1) * (one / *it);
    if (sigma != sigma_) {
      std::cout << sigma << " " << sigma_ << std::endl;
      throw std::runtime_error("No solution found.");
    }
  }

  return gf_polynomial(field, { sigma });
}

gf_polynomial
rs::error_values_naive(const std::vector<gf_element> &syndromes,
                       const std::vector<gf_element> &zeroes) const {
  /* Let s_j be the syndrome values.
   * Let x_i be the error positions.
   * Let y_i be the error values.
   *
   * To find the y_i's solve the linear equation system:
   * 
   * s_j = Σ_{i=1}^{v} y_i * x_i^j
   */

  std::cout << std::endl << "Calculating error values naively" << std::endl;
  const unsigned v = zeroes.size();
  const gf_polynomial factor(
      field, std::vector<gf_element>(std::rbegin(zeroes), std::rend(zeroes)));

  gf_polynomial row(factor);
  gf_matrix system(field);

  
  for(unsigned i = 0; i < v; i++) {
    gf_polynomial tmp(row);
    tmp.push_front(syndromes.at(i));
    system.push_back(tmp);
    // row *= factor; element wise mult.
    std::transform(
        std::begin(row), std::end(row), std::begin(factor), std::begin(row),
        [](const gf_element &lhs, const gf_element &rhs) { return rhs * lhs; });
  }

  std::cout << system << std::endl;

  // TODO Check if the determinant is non-zero.
  std::cout << "Reduced row echelon form of the linear equation system:"
            << std::endl;

  std::cout << system.reduced_echelon_form() << std::endl;

  std::cout << "Solution vector of the linear equation system:" << std::endl;
  auto solution = system.solution();
  std::cout << solution << std::endl << std::endl;

  gf_polynomial e(field, std::vector<int>(n, 0));

  for(int i = 0; i < v; i++) {
    unsigned position = field.polynomial_to_power(zeroes.at(i));
    e[position] = solution[i];
  }

  std::cout << "Error vector e: " << e << std::endl;

  return e;
}

gf_polynomial rs::correct_pzg(const gf_polynomial &b) const {
  std::vector<gf_element> syndromes;

  std::cout << "Syndromes: ";
  for (unsigned power = mu; power < mu + dmin - 1; power++) {
    syndromes.push_back(b(field.power_to_polynomial(power)));
    std::cout << syndromes.back() << " ";
  }
  std::cout << std::endl;

  /* check if all syndrome values are zero - error free code word */
  if (std::all_of(std::begin(syndromes), std::end(syndromes),
                  [&](const auto &e) { return e == field.zero(); })) {
    return b;
  }

  auto sigma = pzg(syndromes);
  sigma.push_back(gf_element(field, 1));
  std::cout << "Σ(x) = ";
  for (size_t i = 0; i < sigma.degree() + 1; i++) {
    std::cout << "x^" << i << " " << sigma[i];
    if (i < sigma.degree())
      std::cout << " + ";
  }
  std::cout << std::endl;
  
  auto lambda = berlekamp_massey(field, syndromes).reverse();
  std::cout << "Λ(x) = ";
  for (size_t i = 0; i < lambda.degree() + 1; i++) {
    std::cout << "x^" << i << " " << lambda[i];
    if (i < lambda.degree())
      std::cout << " + ";
  }
  std::cout << std::endl;

  auto lambda_euklid = ::euklid(field, fk, syndromes);

  if(lambda != sigma) {
    std::cout << "Mistmatch" << std::endl;
    std::cout << "Λ(x) = " << lambda << std::endl;
    std::cout << "Σ(x) = " << sigma << std::endl;
  }

  auto zeroes = sigma.zeroes();
  if (zeroes.size() != sigma.degree()) {
    std::ostringstream os;
    os << "Σ(x) has to have " << sigma.degree() << " zeroes, but it has "
       << zeroes.size() << "." << std::endl;
    throw std::logic_error(os.str());
  }

  std::cout << "Zeroes in Σ(x): ";
  for (const auto &zero : zeroes)
    std::cout << zero << " ";
  std::cout << std::endl;

  auto e = error_values_naive(syndromes, zeroes);
  
  return e + b;
}

gf_polynomial
rs::correct_erasures(const gf_polynomial &b,
                     const std::vector<gf_element> &erasures) const {
  gf_polynomial b_(b);

  /* initialize erasures, in case it was not already done */
  for (const auto &erasure : erasures) {
    b_.at(field.polynomial_to_power(erasure)) = field.zero();
  }

  std::cout << "Erased vector: " << b_ << std::endl;

  std::vector<gf_element> syndromes;

  std::cout << "Syndromes: ";
  for (unsigned power = mu; power < mu + dmin - 1; power++) {
    syndromes.push_back(b(field.power_to_polynomial(power)));
    std::cout << syndromes.back() << " ";
  }
  std::cout << std::endl;

  auto e = error_values_naive(syndromes, erasures);

  return e + b;
}

gf_polynomial rs::correct_bm(const gf_polynomial &b,
                             const std::vector<gf_element> &erasures) const {
  std::vector<gf_element> syndromes;

  std::cout << "Syndromes: ";
  for (unsigned power = mu; power < mu + dmin - 1; power++) {
    syndromes.push_back(b(field.power_to_polynomial(power)));
    std::cout << syndromes.back() << " ";
  }
  std::cout << std::endl;

  /* check if all syndrome values are zero - error free code word */
  const bool syndromes_zero =
      std::all_of(std::begin(syndromes), std::end(syndromes),
                  [&](const auto &e) { return e == field.zero(); });
  if (erasures.size() == 0 && syndromes_zero)
    return b;

  auto lambda = ::berlekamp_massey(field, syndromes, erasures).reverse();
  std::cout << "Λ(x) = ";
  for (size_t i = 0; i < lambda.degree() + 1; i++) {
    std::cout << "x^" << i << " " << lambda[i];
    if (i < lambda.degree())
      std::cout << " + ";
  }
  std::cout << std::endl;

  auto lambda_euklid = ::euklid(field, fk, syndromes, erasures).reverse();

  if (lambda != lambda_euklid) {
    std::cout << "Mistmatch" << std::endl;
    std::cout << "Λ(x)_bm = " << lambda << std::endl;
    std::cout << "Λ(x)_eu = " << lambda_euklid << std::endl;
  }

  auto zeroes = lambda.zeroes();
  if (zeroes.size() != lambda.degree()) {
    std::ostringstream os;
    os << "Σ(x) has to have " << lambda.degree() << " zeroes, but it has "
       << zeroes.size() << "." << std::endl;
    throw std::logic_error(os.str());
  }

  std::cout << "Zeroes in Σ(x): ";
  for (const auto &zero : zeroes)
    std::cout << zero << " ";
  std::cout << std::endl;

  auto e = error_values_naive(syndromes, zeroes);

  return e + b;
}
