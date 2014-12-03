#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>

#include "codes/cyclic.h"
#include "math/linear_equation_system.h"

namespace cyclic {
template <unsigned q, typename Capability,
          typename Sigma = peterson_gorenstein_zierler_tag,
          unsigned N = (1 << q) - 1, typename Coding = division_tag>
class primitive_bch : public cyclic<q, Capability, Sigma, N, Coding> {

  using Base = cyclic<q, Capability, Sigma, N, Coding>;

public:
  using Element = typename Base::Element;
  using Polynomial = typename Base::Polynomial;

private:
  static Polynomial g() {
    Polynomial g({ Element(1) });

    std::vector<unsigned> powers;
    for (unsigned p = 1; p < 2 * static_cast<unsigned>(Base::t); p += 2) {
      powers.push_back(p);
    }

    for (const auto &power : powers) {
      auto roots = minimal_polynomial_roots(power);
      Polynomial m{ Element(1) };
      for (const auto &root : roots) {
        /* build minimal polynomial */
        m *= Polynomial{ Element::from_power(root), Element(1) };
      }
      g = g.lcm(m);
    }
    return g;
  }

  static std::vector<Element> syndromes() {
    std::vector<Element> syndromes;
    Polynomial g_ = g();

    for (unsigned power = 1; power < 2 * Base::t + 1; power++)
      syndromes.push_back(Element::from_power(power));
    return syndromes;
  }

  /*      r - 1
   *       __
   * m_i = || (x - a^{2^j*i})
   *      j = 0
   */
  static std::vector<unsigned> minimal_polynomial_roots(const unsigned &i) {
    std::vector<unsigned> roots{ i };

    for (unsigned r = 1; r < q; r++) {
      auto power = ((1 << r) * i) % Base::n;
      if (power == i)
        break;
      roots.push_back(power);
    }

    auto power = ((1 << q) * i) % Base::n;
    if (power != i) {
      throw std::runtime_error("Cycle seems not be finished after q elements.");
    }

    return roots;
  }

  static std::vector<Element> error_values(const std::vector<Element> &,
                                           const std::vector<Element> &zeroes) {
    return std::vector<Element>(zeroes.size(), Element(1));
  }

  template <
      typename Return_type = typename Base::Element::storage_type,
      typename InputSequence, typename Tag,
      typename std::enable_if<
          !std::is_same<Tag, peterson_gorenstein_zierler_tag>::value>::type * =
          nullptr>
  std::vector<Return_type> correct(const InputSequence &b,
                                   const std::vector<unsigned> &erasures,
                                   Tag) const {
    return Base::template correct<Return_type>(b, erasures);
  }

  template <typename Return_type = typename Base::Element::storage_type,
            typename InputSequence>
  std::vector<Return_type> correct(const InputSequence &b,
                                   const std::vector<unsigned> &erasures,
                                   peterson_gorenstein_zierler_tag) const {
    if (erasures.empty()) {
      return Base::template correct<Return_type>(b, erasures);
    } else {
      if (erasures.size() > 2 * Base::t)
        throw decoding_failure(
            "Number of erasures exceed error correction capability.");

      auto tmp(b);

      std::vector<std::pair<Polynomial, size_t> > results;

      try {
        for (const auto &erasure : erasures)
          tmp.at(erasure) = typename InputSequence::value_type(0);
        results.push_back(Base::template correct_<Return_type>(
            tmp, std::vector<unsigned>{}, Sigma{}));
      }
      catch (const decoding_failure &e) {
        std::cout << e.what() << std::endl;
      }

      try {
        for (const auto &erasure : erasures)
          tmp.at(erasure) = typename InputSequence::value_type(1);
        results.push_back(Base::template correct_<Return_type>(
            tmp, std::vector<unsigned>{}, Sigma{}));
      }
      catch (const decoding_failure &e) {
        std::cout << e.what() << std::endl;
      }

      if (results.empty()) {
        throw decoding_failure("Erasure decoding failed.");
      }

      /* order by number of errors, ascending */
      std::sort(std::begin(results), std::end(results),
                [](const auto &lhs,
                   const auto &rhs) { return lhs.second < rhs.second; });

      const auto &poly = results.front().first;
      std::vector<Return_type> r;
      r.reserve(b.size());
      std::transform(std::cbegin(poly), std::cend(poly), std::back_inserter(r),
                     [](const auto &e) { return Return_type(e); });
      return r;
    }
  }

public:
  primitive_bch() : Base(g(), syndromes(), &error_values) {}

  template <typename Return_type = typename Base::Element::storage_type,
            typename InputSequence>
  std::vector<Return_type> correct(const InputSequence &b,
                                   const std::vector<unsigned> &erasures =
                                       std::vector<unsigned>()) const {
    return correct<Return_type>(b, erasures, Sigma{});
  }
};
}
