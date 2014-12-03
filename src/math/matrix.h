#pragma once

#include <vector>
#include <type_traits>
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
  using std::vector<std::vector<T> >::vector;
  matrix(const std::vector<T> &v)
      : std::vector<std::vector<T> >(1, v), cols(v.size()) {}
  matrix(const size_t rows, const size_t cols_)
      : std::vector<std::vector<T> >(rows, std::vector<T>(cols_, T())),
        cols(cols_) {}

  size_t columns() const { return cols; }
  size_t rows() const { return this->size(); }

  template <typename A, typename R = typename std::common_type<T, A>::type>
  std::vector<R> operator*(const std::vector<A> &rhs) const {
    check_dimension(rhs);

    std::vector<R> result;
    for (const auto &row : *this)
      result.push_back(std::inner_product(std::cbegin(row), std::cend(row),
                                          std::cbegin(rhs), R(0)));

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

