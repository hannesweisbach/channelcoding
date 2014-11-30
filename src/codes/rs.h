#pragma once

#include "codes/cyclic.h"

namespace cyclic {
template <unsigned q, typename Capability,
          typename Sigma = peterson_gorenstein_zierler_tag,
          unsigned N = (1 << q) - 1, typename Coding = division_tag,
          unsigned mu = 1, unsigned step = 1>
class rs : public cyclic<q, Capability, Sigma, N, Coding, naive_tag> {
  using Base = cyclic<q, Capability, Sigma, N, Coding, naive_tag>;

public:
  using Element = typename Base::Element;
  using Polynomial = typename Base::Polynomial;

private:
  static Polynomial g() {
    Polynomial g({ Element(1) });

    for (unsigned i = 0; i < 2 * Base::t; ++i) {
      auto power = mu + i * step;
      Element root = Element::from_power(power);
      g *= Polynomial({ root, Element(1) });
    }

    return g;
  }

  static std::vector<Element> syndromes() {
    std::vector<Element> syndromes;

    for (unsigned i = 0; i < 2 * Base::t; ++i) {
      auto power = mu + i * step;
      Element root = Element::from_power(power);
      syndromes.push_back(root);
    }
    return syndromes;
  }

  std::vector<Element> error_values(const std::vector<Element> &syndromes,
                                    const std::vector<Element> &zeroes,
                                    naive_tag) const {
    /* Let s_j be the syndrome values.
     * Let x_i be the error positions.
     * Let y_i be the error values.
     *
     * To find the y_i's solve the linear equation system:
     *
     * s_j = Σ_{i=1}^{v} y_i * x_i^j
     */

    std::cout << std::endl << "Calculating error values naively" << std::endl;
    const size_t v = zeroes.size();
    const Polynomial factor(
        std::vector<Element>(std::rbegin(zeroes), std::rend(zeroes)));

    Polynomial row(factor);
    math::linear_equation_system<Polynomial> system;

    for (size_t i = 0; i < v; i++) {
      Polynomial tmp({ syndromes.at(i) });
      std::copy(std::cbegin(row), std::cend(row), std::back_inserter(tmp));
      system.push_back(tmp);
      // row *= factor; element wise mult.
      std::transform(
          std::begin(row), std::end(row), std::begin(factor), std::begin(row),
          [](const Element &lhs, const Element &rhs) { return rhs * lhs; });
    }

    std::cout << system << std::endl;

    std::cout << "Solution vector of the linear equation system:" << std::endl;
    auto solution = system.solution();
    std::cout << solution << std::endl << std::endl;

    return solution;
  }
  
  std::vector<Element> error_values(const std::vector<Element> &syndromes,
                                    const std::vector<Element> &zeroes,
                                    forney_tag) const {
    return error_values(syndromes, zeroes, naive_tag());
  }


public:
  rs() : Base(g(), syndromes()) {}

  /* TODO: typedef codeword type? */
  /* InputSequence concepts:
   * size()
   * convertible to bool
   */

};
}
