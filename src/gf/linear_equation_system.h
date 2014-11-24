#pragma once

#include "polynomial.h"

template <typename Polytype>
class linear_equation_system : public std::vector<Polytype> {
  std::vector<Polytype> rows;

  linear_equation_system reduced_echelon_form() const {
    std::vector<Polytype> nrows(*this);

    /* make sure rows are sorted with the left-most elements at the top */
    std::sort(std::begin(nrows), std::end(nrows),
              [](const Polytype &lhs,
                 const Polytype &rhs) { return lhs.degree() > rhs.degree(); });

    for (auto first_row = nrows.begin(); first_row != nrows.end();
         ++first_row) {
      (*first_row) *= first_row->at(first_row->degree()).inverse();
      for (auto next_rows = first_row + 1; next_rows != nrows.end();
           ++next_rows) {
        /* only add if not already zero */
        if (next_rows->at(first_row->degree()))
          (*next_rows) += (*first_row) * next_rows->at(next_rows->degree());
      }
      // std::cout << matrix(nrows) << std::endl;
    }

    for (auto modify = nrows.rbegin() + 1; modify != nrows.rend(); ++modify) {
      for (auto row = nrows.crbegin(); row != modify; ++row) {
        const size_t index = std::distance(std::crbegin(nrows), row) + 1;
        auto factor = modify->at(index);
        (*modify) += (*row) * factor;
      }
    }

    return linear_equation_system(nrows);
  }

public:
  using std::vector<Polytype>::vector;

  linear_equation_system() = default;
  linear_equation_system(const std::vector<Polytype> &v)
      : std::vector<Polytype>::vector(v) {}
  linear_equation_system(std::vector<Polytype> &&v)
      : std::vector<Polytype>::vector(std::move(v)) {}

  Polytype solution() const {
    linear_equation_system copy(reduced_echelon_form());

    /*
     * Matrix:
     * x = 0 0 1
     * y = 0 1 0
     * z = 1 0 0
     *     ^ back()[1].
     * Must be non-zero, otherwise the system is overdetermined.
     * Since it is in row-echelon form, it must be one.
     */
    if (copy.back()[1] != typename Polytype::element_type(1))
      throw std::runtime_error("Linear equation system not solvable");

    Polytype s;
    for (const auto &row : copy) {
      s.push_back(row[0]);
    }

    return s;
  }
};

