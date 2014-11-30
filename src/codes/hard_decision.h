#pragma once

#include <vector>
#include <stdexcept>
#include <iterator>
#include <numeric>
#include <string>

#include "codes/codes.h"
#include "gf/linear_equation_system.h"
#include "protocol.h"

namespace cyclic {

struct peterson_gorenstein_zierler_tag : hard_decision_tag {
  static std::string to_string() { return "PGZ"; }
};
struct berlekamp_massey_tag : hard_decision_tag {
  static std::string to_string() { return "BM"; }
};
/* AKA Sugiyama algorithm */
struct euklid_tag : hard_decision_tag {
  static std::string to_string() { return "EUKLID"; }
};
/* TODO:
 * L. R. Welch, E. Berlekamp, Error correction for algebraic block codes, US
 *                            Patent, Number 4,633,470, 1986.
 * G. Forney, Generalized minimum distance decoding, Information Theory, IEEE
 *            Transactions on, vol.12, no.2, pp. 125-131, 1966.
 * J. Jiang, K. R. Narayanan, Iterative Soft Decoding of Reed-Solomon Codes,
 *                            IEEE Communications Letters, vol.8, no.4,
 *                            pp.244-246, 2004.
 * inverse-free BMAs
 */

#if 0
std::vector<gf_element> chien(const gf &field, gf_polynomial polynomial) {
  std::vector<gf_element> zeroes;
  std::vector<gf_element> factors(polynomial.size(), field.zero());
  std::iota(std::begin(factors), std::end(factors), field.one());

  for (int i = 0; i < field.size(); i++) {
    auto result = std::accumulate(std::begin(polynomial), std::end(polynomial),
                                  field.zero());
    if (result == field.zero()) {
      zeroes.push_back(field.power_to_polynomial(i));
    }

    /* multiply by alpha, alpha^2, ... */
    std::transform(
        std::cbegin(polynomial), std::cend(polynomial), std::cbegin(factors),
        std::begin(polynomial),
        [](const gf_element &lhs, const gf_element &rhs) { return lhs * rhs; });
  }

  return zeroes;
}
#endif


template <typename Polynomial,
          typename Element = typename Polynomial::coefficient_type>
Polynomial error_locator_polynomial(const std::vector<Element> &syndromes,
                                    const std::vector<unsigned> &erasures,
                                    peterson_gorenstein_zierler_tag) {
  if (!erasures.empty())
    throw std::runtime_error(
        "The PGZ-Algorithm does not support erasure decoding");

  /* TODO: this implementation suffers from decoder malfunction. Fix it.
   * http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=333881
   */
  for (ssize_t v = static_cast<ssize_t>(syndromes.size() / 2); v; --v) {
    linear_equation_system<Polynomial> eq_system;
    for (auto it = std::cbegin(syndromes); it != std::cbegin(syndromes) + v;
         ++it) {
      Polynomial poly;
      std::reverse_copy(it, it + v + 1, std::back_inserter(poly));
      eq_system.push_back(poly);
    }

    // TODO Check if the determinant is non-zero.
    try {
      auto solution = eq_system.solution();
      solution.push_back(Element(1));
      return solution;
    }
    catch (const std::exception &) {
    }
  }

  /* note: this was supposed to solve decoder malfunction, but it is bullshit.
   * if v < t, then multiple equation systems can be built. each must be
   * solvable and yield the same result. otherwise, decoding failure must be
   * declared.
   * This is not an immediate problem anymore, since the corrected codeword's
   * syndrome is checked again.
   */
  auto sigma = syndromes.front();
  std::cout << sigma << std::endl;
  for (auto it = std::begin(syndromes); it != std::end(syndromes) - 1; ++it) {
    if (*it == Element(0))
      throw decoding_failure("Cannot invert syndrome, it is 0.");
    std::cout << *it << " " << Element(1) / *it << std::endl;
    auto sigma_ = *(it + 1) * (Element(1) / *it);
    std::cout << *(it + 1) << " " << sigma_ << std::endl;
    if (sigma != sigma_) {
      // std::cout << sigma << " " << sigma_ << std::endl;
      throw decoding_failure("No solution found.");
    }
  }

  return Polynomial({ sigma, Element(1) });
}

template <typename Polynomial,
          typename Element = typename Polynomial::coefficient_type>
Polynomial error_locator_polynomial(const std::vector<Element> &syndromes,
                                    const std::vector<unsigned> &erasures,
                                    berlekamp_massey_tag) {
  const auto fk = syndromes.size() / 2;
  const auto rho = erasures.size();
  Polynomial lambda({ Element(1) });
  Polynomial b;
  auto l = erasures.size();

  /* position of the erasure is the power of the element */
  for (const auto &erasure : erasures)
    lambda *= Polynomial({ Element(1), Element::from_power(erasure) });

  b = lambda;

  for (auto i = rho; i < 2 * fk; i++) {
    /* b = b * x; */
    b *= Polynomial({ Element(0), Element(1) });
    /* d = si + \sigma_j=1^l lambda_j * s_i-j; */
    ssize_t end_offset = static_cast<ssize_t>(l + 1);
    ssize_t start_offset = static_cast<ssize_t>(i);
    const auto delta = std::inner_product(
        std::cbegin(lambda) + 1, std::cbegin(lambda) + end_offset,
        std::crend(syndromes) - start_offset, syndromes.at(i));

    if (delta) {
      auto t = lambda + b * delta;
      if (2 * l <= i + rho) {
        b = lambda * delta.inverse();
        l = i + rho - l + 1;
      }
      lambda = t;
    }
    protocol_bm(i, delta, lambda, l, b, std::false_type());
  }

  return lambda.reverse();
}

template <typename Polynomial,
          typename Element = typename Polynomial::coefficient_type>
Polynomial error_locator_polynomial(const std::vector<Element> &syndromes,
                                    const std::vector<unsigned> &erasures,
                                    euklid_tag) {
  const auto fk = syndromes.size() / 2;
  const auto dmin = 2 * fk + 1;
  const auto max = static_cast<ssize_t>((2 * fk + erasures.size()) / 2);

  Polynomial u({ Element(1) });
  std::vector<Polynomial> w;
  std::vector<Polynomial> r;
  Polynomial s(syndromes);

  for (const auto &erasure : erasures)
    u *= Polynomial({ Element(1), Element::from_power(erasure) });

  r.push_back(s * u);
  r.push_back(Polynomial(dmin - 1, Element(0)));
  r.back().push_back(Element(1));

  w.push_back(u);
  w.push_back(Polynomial({ Element(0) }));

  for (size_t i = 1; r.back().degree() >= max; i++) {
    auto next = r.at(i - 1) % r.at(i);
    auto q = r.at(i - 1) / r.at(i);
    w.push_back(w.at(i - 1) + q * w.at(i));
    r.push_back(next);
  }

  /* I wish I had static_if ... */
  protocol_euklid(r, w, std::false_type());

  if (w.back().at(0) == Element(0))
    throw decoding_failure("Cannot invert last element");

  auto lambda = w.back() * w.back().at(0).inverse();
  return lambda.reverse();
}
}
