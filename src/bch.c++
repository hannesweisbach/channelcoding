#include <stdexcept>
#include <iterator>
#include <utility>
#include <typeinfo>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <tuple>

#include <iomanip>
#include <cstdlib>
#include <clocale>

#include "util.h"
#include "bch.h"
#include "cyclic.h"
#include "center.h"

std::vector<gf_element> chien(const gf &field,
                              gf_polynomial polynomial) {
  std::vector<gf_element> zeroes;
  std::vector<gf_element> factors(polynomial.size(), field.zero());
  std::iota(std::begin(factors), std::end(factors), field.one());

  for(int i = 0; i < field.size(); i++) {
    auto result = std::accumulate(std::begin(polynomial), std::end(polynomial),
                                  field.zero());
    if(result == field.zero()) {
      zeroes.push_back(field.power_to_polynomial(i));
    }

    /* multiply by alpha, alpha^2, ... */
    std::transform(
        std::cbegin(polynomial), std::cend(polynomial), std::cbegin(factors)
        , std::begin(polynomial),
        [](const gf_element &lhs, const gf_element &rhs) { return lhs * rhs; });
  }

  return zeroes;
}

static unsigned consecutive_powers(std::vector<unsigned> powers) {
  std::sort(std::begin(powers), std::end(powers));
  unsigned l = 1;
  
  for(auto it = std::cbegin(powers); it != std::cend(powers) - 1; ++it) {
    if (((*it) + 1) == *(it + 1))
      l++;
    else
      //TODO start new count
      break;
  }

  return l;
}

static unsigned int next_missing(const std::vector<unsigned int> powers) {
  //TODO: assert is sorted.
  for (unsigned int i = 0; i < powers.size(); i++) {
    if (powers[i] != i)
      return i;
  }

  return powers.size();
}

std::vector<std::vector<unsigned int> > bch_cycles(unsigned int n) {
  std::vector<unsigned int> powers(1, 0);
  std::vector<std::vector<unsigned int> > cycles(
      1, std::vector<unsigned int>({ 0 }));

  for (; powers.size() < n;) {
    std::sort(std::begin(powers), std::end(powers));
    unsigned int i = next_missing(powers);
    std::vector<unsigned int> cycle(1, i);
    powers.push_back(i);

    for(i = (i * 2) % n; i != cycle[0]; i = (i*2) %n) {
      cycle.push_back(i);
      powers.push_back(i);
    }

    std::sort(std::begin(cycle), std::end(cycle));
    cycles.push_back(cycle);
  }

  return cycles;
}

bch::bch(unsigned q, uint64_t modular_polynomial, unsigned d_e, enum type type)
    : bch(q, modular_polynomial, d_e, (1 << q) - 1) {
  this->type = type;
}

bch::bch(unsigned q, uint64_t modular_polynomial, unsigned d_e, unsigned n)
    : field(q, modular_polynomial), n(n), l(0), k(0), dmin(0),
      type(type::shortened), generator(field, { 1 }), f(field, { 1 }),
      h(field), control_matrix(0, 0) {
  if (d_e > n)
    throw std::logic_error("d_e <= n required");
  // k = degree(g(x)) = sum_µ^µ+d-2 deg(m_i(x));
  // l = n - k;
  // dmin = # of consecutive zeroes / powers of alpha in g(x) + 1
  std::cout << "BCH code of length n = 2^q - 1 = " << n << std::endl;
  std::cout << "k = deg(g(x)) = Σ_µ^{µ+d-2} deg(m_i(x)) = " << std::endl;
  std::cout << "l = n - k = " << n << " - " << k << " = " << l << std::endl;
  auto cycles = bch_cycles(n);
  for (const auto &cycle : cycles) {
    std::cout << "{ ";
    std::copy(std::begin(cycle), std::end(cycle),
              std::ostream_iterator<unsigned int>(std::cout, ", "));
    // for(const auto &power : cycle)
    //  std::cout << power  << ", ";
    std::cout << "}; degree: " << cycle.size() << std::endl;
  }

  std::vector<std::tuple<int, int, int, int>> codes;
  for (auto it = std::begin(cycles); it != std::end(cycles); ++it) {
    mu = it->front();
    std::vector<unsigned> powers;
    k = 0;
    int number_of_cycles = 0;
    for (auto cycle = it; cycle != std::end(cycles); ++cycle) {
      powers.insert(std::end(powers), std::begin(*cycle), std::end(*cycle));
      dmin = consecutive_powers(powers) + 1;

      if (dmin > n)
        break;

      number_of_cycles++;
      k += cycle->size();

      auto code =
          std::find_if(std::begin(codes), std::end(codes),
                       [&](const auto &c) { return std::get<0>(c) == dmin; });

      if (code != std::end(codes)) {
        /* code found */
        if (k < std::get<1>(*code)) {
          /* overwrite with better code */
          *code = std::make_tuple(dmin, k, mu, number_of_cycles);
        }
        break;
      }

      codes.emplace_back(dmin, k, mu, number_of_cycles);
    }
  }

  std::vector<std::tuple<int, int, int, int>> matching_codes;
  for (const auto &code : codes) {
    std::tie(dmin, k, mu, std::ignore) = code;
    
    if (dmin >= d_e)
      matching_codes.push_back(code);
    
    l = n - k;
    std::cout << "(" << n << ", " << l << ", " << dmin
              << ") BCH code found with µ = " << mu << std::endl;
  }

  if (matching_codes.empty())
    throw std::runtime_error("No codes found to fulfill requirements.");

  std::cout << std::endl << "Codes with dmin >= d_e:" << std::endl;

  for(const auto&code : matching_codes) {
    std::tie(dmin, k, mu, std::ignore) = code;
    l = n - k;
    std::cout << "(" << n << ", " << l << ", " << dmin
              << ") BCH code found with µ = " << mu << std::endl;

  }

  std::sort(std::begin(matching_codes), std::end(matching_codes),
            [](const auto &lhs, const auto &rhs) {
    return std::get<1>(lhs) < std::get<1>(rhs);
  });

  int required_cycles;
  const auto & code = matching_codes.front();
  std::tie(dmin, k, mu, required_cycles) = code;

  l = n - k;
  std::cout << "Selected code: (" << n << ", " << l << ", " << dmin
            << ") BCH, µ = " << mu << " to " << required_cycles << std::endl;
  std::cout << std::endl;

  for (auto cycle = std::begin(cycles) + mu;
       cycle != std::begin(cycles) + mu + required_cycles; ++cycle) {
    for (const auto &power : *cycle) {
      generator *= gf_polynomial(field, { field.power_to_polynomial(power),
                                          field.power_to_polynomial(0) });
    }
  }
  std::cout << "Generator polynomial g(x) = " << generator << std::endl;

  //f(x) = x^n - 1; 
  const gf_polynomial x(field, std::vector<int>({ 0, 1 }));
  for(unsigned i = 0; i < n; ++i)
    f *= x;
  f += gf_polynomial(field, { 1 });
  std::cout << "f(x) = " << f << std::endl;

  h = f / generator;
  std::cout << "h(x) = " << h << std::endl;

  fk = (dmin - 1) / 2;

  /* control matrix H:
   *     ,              ,
   *     | 0   …   f(x) |
   *     |       …      |
   * H = | 0 … f(x) … 0 |
   *     |       …      |
   *     |f(x)      …   |
   *     `              ´
   */

  auto h_tmp(h);
  h_tmp.reverse();
  
  for(unsigned i = 0; i < k; ++i) {
    std::vector<int> tmp;
    for(const auto & coeff : h_tmp)
      tmp.push_back(coeff ? 1 : 0);
    // add current f
    control_matrix.push_back(tmp);
    // cyclic rotation
    h_tmp *= x;
  }
  
  control_matrix.resize(k, n);
  //std::reverse(std::begin(control_matrix), std::end(control_matrix));
}

const matrix<int> &bch::H() const { return control_matrix; }

std::vector<int> bch::encode_div(const std::vector<int> &a) const {
  if (a.size() > l) {
    std::ostringstream os;
    os << "Input word (" << a.size() << ") longer than l = " << l;
    throw std::runtime_error(os.str());
  }

  gf_polynomial b(field, a);
  const gf_polynomial x(field, std::vector<int>({ 0, 1 }));

  std::generate_n(std::back_inserter(b), l - a.size(),
                  [&]() { return field.zero(); });
  for (unsigned i = 0; i < k; i++)
    b *= x;

  auto r = b % generator;
  b += r;

  std::cout << "Code word: " << b << std::endl;
  
  std::vector<int> result;
  for (const auto &e : b)
    result.push_back((e) ? 1 : 0);

  return result;
}

std::vector<gf_element> bch::syndromes(const gf_polynomial &b) const {
  std::vector<gf_element> syndrome;
  const unsigned fk = (dmin - 1) / 2;
  unsigned t = fk * 2;

  if(t != mu + dmin - 2) {
    std::ostringstream os;
    os << "t != µ + dmin - 2. t = " << t << ", µ = " << mu
       << ", dmin = " << dmin << ", µ + dmin - 2 = " << mu + dmin - 2
       << std::endl;
    throw std::runtime_error(os.str());
  }

  //std::cout << "Syndrome values: ";
  for (unsigned i = mu; i <= mu + dmin - 2; i++) {
    syndrome.push_back(b(field.power_to_polynomial(i)));
    //std::cout << syndrome.back() << " ";
  }

  //std::cout << std::endl;

  return syndrome;
}

gf_polynomial bch::pzg(std::vector<gf_element> syndromes) const {
  int s0 = fk % 2;

  if (mu == 0) {
    s0 = syndromes.front();
    if (!((s0 == field.zero()) || (s0 == field.one()))) {
      std::ostringstream os;
      os << "s0 has to be either 0 or 1. It was " << s0 << std::endl;
      throw std::logic_error(os.str());
    }
    syndromes.erase(std::begin(syndromes));
  }

  const int x = ((fk % 2) ^ s0);
  for (unsigned v = fk - x; v > x; v = v - (1 + x)) {
    gf_matrix eq_system(field);
    for (auto it = std::cbegin(syndromes); it != std::cbegin(syndromes) + v;
         ++it) {
      gf_polynomial poly(field);
      for (auto cp_it = it; cp_it != it + v + 1; ++cp_it)
        poly.push_front(*cp_it);
      eq_system.push_back(poly);
    }

    //std::cout << "Linear equation system to be solved: " << std::endl;
    //std::cout << eq_system << std::endl;

    // TODO Check if the determinant is non-zero.
    try {
#if 0
      std::cout << "Reduced row echelon form of the linear equation system:"
                << std::endl;

      std::cout << eq_system.reduced_echelon_form() << std::endl;

      std::cout << "Solution vector of the linear equation system:"
                << std::endl;
#endif
      auto solution = eq_system.solution();
      //std::cout << solution << std::endl << std::endl;

      return solution;
    }
    catch (const std::exception &e) {
#if 0
      std::cout << "Error for v = " << v << ": " << typeid(e).name() << " "
                << e.what() << std::endl;
#endif
    }
  }

  /* if v = 1, check all equations for the same solution for sigma_1 */
  auto sigma = syndromes.front();
  for (auto it = std::begin(syndromes); it != std::end(syndromes) - 1; ++it) {
    auto sigma_ = *(it + 1) * (field.one() / *it);
    if (sigma != sigma_) {
      std::cout << sigma << " " << sigma_ << std::endl;
      throw std::runtime_error("No solution found.");
    }
  }

  return gf_polynomial(field, { sigma });
}

std::vector<int> bch::correct_peterson(const std::vector<int> &b_) const {

  gf_polynomial b(field, b_);

  std::vector<gf_element> syndrome = syndromes(b);

  /* check if all syndrome values are zero - error free code word */
  auto it = std::find_if(std::begin(syndrome), std::end(syndrome),
                         [&](const auto &e) { return e != field.zero(); });
  if (it == std::end(syndrome))
    return b_;


  auto solution = pzg(std::move(syndrome));
  solution.push_back(gf_element(field, 1));
#if 0
  std::cout << "Σ(x) = ";
  for (size_t i = 0; i < solution.degree() + 1; i++) {
    std::cout << "x^" << i << " " << solution[i];
    if (i < solution.degree())
      std::cout << " + ";
  }
  std::cout << std::endl;
#endif
  
  //chien(field, solution);
  auto zeroes = solution.zeroes();
  if (zeroes.size() != solution.degree()) {
    std::ostringstream os;
    os << "Σ(x) has to have " << solution.degree() << " zeroes, but it has "
       << zeroes.size() << "." << std::endl;
    throw std::logic_error(os.str());
  }

  std::vector<int> e(n, 0);

#if 0
  std::cout << "Zeroes in Σ(x): ";
  for (const auto &zero : zeroes) {
    std::cout << zero << " ";
    e.at(zero.power()) = 1;
  }
  std::cout << std::endl;

  std::cout << "Error vector     e (x) = ";
  for (const auto &i : e)
    std::cout << i << " ";
  std::cout << std::endl;
#endif
  
  // Corrected vector in e (b_ is const)
  std::transform(std::begin(e), std::end(e), std::begin(b_), std::begin(e),
                 [](const int &lhs, const int &rhs) { return lhs ^ rhs; });
#if 0
  std::cout << "Received vector  b (x) = ";
  for (const auto &i : b_)
    std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "Corrected vector b*(x) = ";
  for (const auto &i : e)
    std::cout << i << " ";
  std::cout << std::endl;
#endif
  return e;

  throw std::runtime_error("No viable solution found. Decoding failure");
}

std::vector<int> bch::decode(const std::vector<int> &b) const {
  if (b.size() == n + 1) {

  } else {
  }
}

gf_polynomial bch::correct_bm(const gf_polynomial &b,
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

  auto lambda = berlekamp_massey(field, syndromes, erasures).reverse();
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

  //auto e = error_values_naive(syndromes, zeroes);

  return b;
}

