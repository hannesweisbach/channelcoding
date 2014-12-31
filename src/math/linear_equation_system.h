#pragma once

namespace math {

template <typename Polytype>
class linear_equation_system {
  using row_type = std::vector<Polytype>;
  using iterator = typename row_type::iterator;
  using const_iterator = typename row_type::const_iterator;
  std::vector<Polytype> rows;

  linear_equation_system reduced_echelon_form() const {
    std::vector<Polytype> nrows(this->rows);

    /* make sure rows are sorted with the left-most elements at the top */
    std::sort(std::begin(nrows), std::end(nrows),
              [](const Polytype &lhs,
                 const Polytype &rhs) { return lhs.degree() > rhs.degree(); });

    /* empty rows */
    if (nrows.front().degree() < 0)
      return linear_equation_system(nrows);

    for (auto first_row = nrows.begin(); first_row != nrows.end();
         ++first_row) {
      size_t offset = static_cast<size_t>(first_row->degree());
      (*first_row) *= first_row->at(offset).inverse();
      for (auto next_rows = first_row + 1; next_rows != nrows.end();
           ++next_rows) {
        /* only add if not already zero */
        if (next_rows->at(offset)) {
          size_t next_offset = static_cast<size_t>(next_rows->degree());
          (*next_rows) += (*first_row) * next_rows->at(next_offset);
        }
      }
      // std::cout << matrix(nrows) << std::endl;
    }

    for (auto modify = nrows.rbegin() + 1; modify != nrows.rend(); ++modify) {
      for (auto row = nrows.crbegin(); row != modify; ++row) {
        const size_t index =
            static_cast<size_t>(std::distance(std::crbegin(nrows), row) + 1);
        auto factor = modify->at(index);
        (*modify) += (*row) * factor;
      }
    }

    return linear_equation_system(std::move(nrows));
  }

public:
  linear_equation_system() = default;
  linear_equation_system(const std::vector<Polytype> &v) : rows(v) {}
  linear_equation_system(std::vector<Polytype> &&v) : rows(std::move(v)) {}
  linear_equation_system(linear_equation_system &&) = default;

  void push_back(const Polytype &row) { rows.push_back(row); }
  void push_back(Polytype &&row) { rows.push_back(std::move(row)); }
  Polytype &back() { return rows.back(); }
  const Polytype &back() const { return rows.back(); }

  iterator begin() noexcept { return rows.begin(); }
  const_iterator begin() const noexcept { return rows.begin(); }
  iterator end() noexcept { return rows.end(); }
  const_iterator end() const noexcept { return rows.end(); }

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
    if (copy.back()[1] != typename Polytype::coefficient_type(1))
      throw std::runtime_error("Linear equation system not solvable");

    Polytype s;
    for (const auto &row : copy) {
      s.push_back(row[0]);
    }

    return s;
  }

  template <class charT, class traits>
  friend std::basic_ostream<charT, traits> &operator<<(
      std::basic_ostream<charT, traits> &os,
      const linear_equation_system &les) {
    for (const auto &row : les) {
      os << row << std::endl;
    }
    return os;
  }
};
}
