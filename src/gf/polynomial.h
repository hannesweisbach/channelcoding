#pragma once

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>

namespace gf {

struct brute_force_tag {};
struct chien_tag {};

template <typename Field, typename Coefficient = typename Field::element_type,
          typename Polynomial>
std::vector<Coefficient> roots(const Polynomial &p, brute_force_tag) {
  /* I need a type for the variable x, it has not necessarily something to do
   * with the coefficient type */
  std::vector<Coefficient> zeroes;

  for (const auto &element : Field{}) {
    if (element && !p(element))
      zeroes.push_back(element);
  }
  return zeroes;
}
}
