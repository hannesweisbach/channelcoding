#pragma once

#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>

template <typename T> class matrix;
template <typename T>
std::ostream &operator<<(std::ostream &, const matrix<T> &);

template <typename T> class matrix : public std::vector<std::vector<T> > {
  size_t cols;

  template <typename A> void check_dimension(const std::vector<A> &rhs) const {
    if (rhs.size() != cols) {
      std::ostringstream os;
      os << "Vector dimension " << rhs.size()
         << " does not match matrix dimension (" << rows() << " x " << columns()
         << ")";
      throw std::runtime_error(os.str());
    }
  }

public:
  matrix(const size_t rows, const size_t cols)
      : std::vector<std::vector<T> >(rows, std::vector<T>(cols, T())),
        cols(cols) {}

  size_t columns() const { return cols; }
  size_t rows() const { return this->size(); }

  void resize(size_t nrows, size_t ncols) {
    /* treat existing rows */
    for (auto &&row : *this) {
      /* erase */
      const auto zero = typename std::vector<std::vector<T> >::size_type();
      row.erase(std::end(row) - std::min(ncols - row.size(), zero),
                std::end(row));
      /* and fill up */
      std::fill_n(std::back_inserter(row), ncols - row.size(), T());
    }
    /* fill up rows to nrows with empty rows */
    std::fill_n(std::back_inserter(*this), nrows - rows(),
                std::vector<T>(ncols, T()));

    cols = ncols;
  }

  template <typename A,
            typename R = decltype(std::declval<T>() * std::declval<A>())>
  std::vector<R> operator*(const std::vector<A> &rhs) const {
    check_dimension(rhs);

    std::vector<R> result;
    for (const auto &row : *this)
      result.push_back(std::inner_product(std::cbegin(row), std::cend(row),
                                          std::cbegin(rhs), T()));

    return result;
  }

  /* wedge product */
  template <typename A,
            typename R = decltype(std::declval<T>() * std::declval<A>())>
  matrix<R> wedge_product(const std::vector<A> &rhs) const {
    check_dimension(rhs);

    matrix<R> result(rows(), columns());

    for (size_t row = 0; row < rows(); row++) {
      for (size_t col = 0; col < columns(); col++) {
        result.at(row).at(col) = this->at(row).at(col) * rhs.at(col);
      }
    }

    return result;
  }

  friend std::ostream &operator<<<>(std::ostream &os, const matrix<T> &matrix);
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const matrix<T> &matrix) {
  for (const auto &row : matrix) {
    for (const auto &e : row)
      os << std::setw(2) << e << " ";
    os << std::endl;
  }
  return os << std::endl;
}

