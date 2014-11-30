#pragma once

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>

namespace gf {

struct brute_force_tag {};
struct chien_tag {};

template <typename Polynomial,
          typename Coefficient = typename Polynomial::coefficient_type>
std::vector<Coefficient> roots(const Polynomial &p, brute_force_tag) {
  typename Coefficient::Field field;
  std::vector<Coefficient> zeroes;

  for (const auto &element : field)
    if (p(element) == Coefficient(0))
      zeroes.push_back(element);
  return zeroes;
}

}
