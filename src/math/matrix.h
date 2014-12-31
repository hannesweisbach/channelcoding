#pragma once

#include <vector>
#include <type_traits>
#include <iomanip>

template <typename T> class matrix;
template <typename T>
std::ostream &operator<<(std::ostream &, const matrix<T> &);

template <typename T> class matrix {
  using row_type = std::vector<T>;
  using rep_type = std::vector<row_type>;

  using iterator = typename rep_type::iterator;
  using const_iterator = typename rep_type::const_iterator;

  rep_type data;
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
  using value_type = row_type;
  using size_type = typename row_type::size_type;

  matrix() = default;
  matrix(const std::vector<T> &v) : data(1, v), cols(v.size()) {}
  matrix(const size_t rows, const size_t cols_)
      : data(rows, std::vector<T>(cols_, T())), cols(cols_) {}
  matrix(const matrix &) = default;
  matrix(matrix &&) = default;

  void push_back(const row_type &e) { data.push_back(e); }
  void push_back(row_type &&e) { data.push_back(std::move(e)); }
  row_type &at(size_type i) { return data.at(i); }
  const row_type &at(size_type i) const { return data.at(i); }

  iterator begin() noexcept { return data.begin(); }
  iterator end() noexcept { return data.end(); }
  const_iterator begin() const noexcept { return data.begin(); }
  const_iterator end() const noexcept { return data.begin(); }
  const_iterator cbegin() const noexcept { return data.cbegin(); }
  const_iterator cend() const noexcept { return data.cend(); }

  size_t columns() const { return cols; }
  size_t rows() const { return data.size(); }

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

