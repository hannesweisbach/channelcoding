#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>

#include "gf/gf.h"
#include "gf/polynomial.h"
#include "codes/cyclic.h"
#include "gf/linear_equation_system.h"

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

  std::vector<Element> error_values(const std::vector<Element> &,
                                    const std::vector<Element> &zeroes,
                                    naive_tag) const override {
    return std::vector<Element>(zeroes.size(), Element(1));
  }

  std::vector<Element> error_values(const std::vector<Element> &s,
                                    const std::vector<Element> &z,
                                    forney_tag) const override {
    return error_values(s, z, naive_tag());
  }

public:
  primitive_bch() : Base(g(), syndromes()) {}
  virtual ~primitive_bch() = default;

  primitive_bch(const primitive_bch &) = default;
  primitive_bch(primitive_bch &&) = default;
  primitive_bch &operator=(const primitive_bch &) = default;
  primitive_bch &operator=(primitive_bch &&) = default;

  template <typename Return_type = typename Base::Galois_Field::storage_t,
            typename InputSequence>
  std::vector<Return_type> correct(const InputSequence &b,
                                   const std::vector<unsigned> &erasures =
                                       std::vector<unsigned>()) const {
    if (!std::is_same<Sigma, peterson_gorenstein_zierler_tag>::value ||
        erasures.empty()) {
      return Base::template correct<Return_type>(b, erasures);
    } else {
      if (erasures.size() > 2 * Base::t)
        throw decoding_failure(
            "Number of erasures exceed error correction capability.");

      auto tmp(b);

      std::vector<Return_type> v0;
      std::vector<Return_type> v1;

      try {
        for (const auto &erasure : erasures)
          tmp.at(erasure) = Return_type(0);
        v0 = Base::template correct<Return_type>(tmp);
      }
      catch (const decoding_failure &e) {
        std::cout << e.what() << std::endl;
      }

      try {
        for (const auto &erasure : erasures)
          tmp.at(erasure) = Return_type(1);
        v1 = Base::template correct<Return_type>(tmp);
      }
      catch (const decoding_failure &e) {
        std::cout << e.what() << std::endl;
      }

      if (v0.empty() && v1.empty()) {
        throw decoding_failure("Erasure decoding failed.");
      } else if (v0.empty()) {
        return v1;
      } else if (v1.empty()) {
        return v0;
      } else {
        auto v0_it = std::cbegin(v0);
        size_t v0_equals =
            std::count_if(std::cbegin(b), std::cend(b),
                          [&](const auto &e) { return e == *v0_it++; });
        auto v1_it = std::cbegin(v1);
        size_t v1_equals =
            std::count_if(std::cbegin(b), std::cend(b),
                          [&](const auto &e) { return e == *v1_it++; });
        std::cout << v0_equals << " " << v1_equals << std::endl;
        return (v0_equals >= v1_equals) ? v0 : v1;
      }
    }
  }
};
}
