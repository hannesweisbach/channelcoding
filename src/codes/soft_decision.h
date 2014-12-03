#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <numeric>
#include <cmath>
#include <tuple>
#include <functional>
#include <type_traits>
#include <ratio>

#include "codes.h"
#include "math/matrix.h"

#ifdef NDEBUG
#define at(x) operator[](x)
#endif

template <unsigned Iterations = 50> struct min_sum_tag : soft_decision_tag {
  static constexpr unsigned iterations = Iterations;
  static std::string to_string() { return "MS"; }
};

namespace detail {
template <typename T> struct is_ratio : std::false_type {
  static constexpr bool value = false;
};

template <std::intmax_t Num, std::intmax_t Denom>
struct is_ratio<std::ratio<Num, Denom> > : std::true_type {
  static constexpr bool value = true;
};
}

template <unsigned Iterations, typename T = std::ratio<1> >
struct normalized_min_sum_tag : soft_decision_tag {
  static_assert(detail::is_ratio<T>::value, "needs to be std::ratio<>.");
  static constexpr unsigned iterations = Iterations;
  static constexpr double alpha = static_cast<double>(T::num) / T::den;
  static std::string to_string() { return "NMS"; }
};

template <unsigned Iterations = 50, typename T = std::ratio<0> >
struct offset_min_sum_tag : soft_decision_tag {
  static_assert(detail::is_ratio<T>::value, "needs to be std::ratio<>.");
  static constexpr unsigned iterations = Iterations;
  static constexpr double beta = static_cast<double>(T::num) / T::den;
  static std::string to_string() { return "OMS"; }
};

template <unsigned Iterations = 50>
struct self_correcting_1_min_sum_tag : soft_decision_tag {
  static constexpr unsigned iterations = Iterations;
  static std::string to_string() { return "SCMS1"; }
};

template <unsigned Iterations = 50>
struct self_correcting_2_min_sum_tag : soft_decision_tag {
  static constexpr unsigned iterations = Iterations;
  static std::string to_string() { return "SCMS2"; }
};

template <unsigned Iterations = 50, typename Alpha = std::ratio<1>,
          typename Beta = std::ratio<1, 10> >
struct normalized_2d_min_sum_tag : soft_decision_tag {
  static_assert(detail::is_ratio<Alpha>::value, "needs to be std::ratio<>.");
  static_assert(detail::is_ratio<Beta>::value, "needs to be std::ratio<>.");
  static constexpr unsigned iterations = Iterations;
  static constexpr double alpha = static_cast<double>(Alpha::num) / Alpha::den;
  static constexpr double beta = static_cast<double>(Beta::num) / Alpha::den;
  static std::string to_string() { return "2DNMS"; }
};

template <typename T> constexpr int signum(const T &val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T, typename U>
inline auto syndrome(const matrix<T> &H, const std::vector<U> &b) {
  const auto syndromes = H * b;
  for (const auto &s : syndromes)
    if (s % 2)
      return false;
  return true;
}

template <typename U, typename R>
std::vector<R> column_sum(const matrix<U> &H, const matrix<R> &r) {
  std::vector<R> col_sums(H.columns(), R(0));
  for (size_t row = 0; row < H.rows(); row++) {
    for (size_t col = 0; col < H.columns(); col++) {
      if (H.at(row).at(col)) {
        col_sums.at(col) += r.at(row).at(col);
      }
    }
  }

  return col_sums;
}

/* check node update */
template <typename U, typename Q, typename R, typename Functor>
void horizontal__(const matrix<U> &H, const matrix<Q> &q, matrix<R> &r,
                  Functor &&fn) {
  const size_t rows = H.rows();
  const size_t cols = H.columns();
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < cols; col++) {
      if (H.at(row).at(col)) {
        int sign = 1;
        Q min = std::numeric_limits<Q>::max();

        for (size_t i = 0; i < cols; i++) {
          if (i != col && H.at(row).at(i)) {
            sign *= signum(q.at(row).at(i));
            min = std::min(min, std::abs(q.at(row).at(i)));
          }
        }
        r.at(row).at(col) = static_cast<R>(sign * fn(min));
      }
    }
  }
}

/* symbol node update */
template <typename U, typename Q, typename R, typename Functor>
void vertical__(const matrix<U> &H, const std::vector<Q> y, const matrix<R> &r,
                matrix<Q> &q, Functor &&fn) {
  std::vector<R> col_sums(column_sum(H, r));

  const size_t rows = H.rows();
  const size_t cols = H.columns();
  for (size_t row = 0; row < rows; row++) {
    for (size_t col = 0; col < cols; col++) {
      if (H.at(row).at(col)) {
        const auto exclusive_colsum = col_sums.at(col) - r.at(row).at(col);
        q.at(row).at(col) = fn(exclusive_colsum, y.at(col), q.at(row).at(col));
      }
    }
  }
}

template <typename U, typename Q, typename R,
          typename Result_t = typename std::common_type<Q, R>::type>
std::vector<Result_t> likelihood(const std::vector<Q> &y, const matrix<U> &H,
                                 const matrix<R> &r) {
  std::vector<Result_t> result;
  result.reserve(y.size());
  std::copy(std::cbegin(y), std::cend(y), std::back_inserter(result));

  for (size_t row = 0; row < H.rows(); row++) {
    for (size_t col = 0; col < H.columns(); col++) {
      if (H.at(row).at(col)) {
        result.at(col) += r.at(row).at(col);
      }
    }
  }

  return result;
}

template <typename U = unsigned, typename T>
std::vector<U> hard_decision(const std::vector<T> &L) {
  static_assert(!std::is_signed<U>::value, "U has to be unsigned");
  std::vector<U> b_hard;
  b_hard.reserve(L.size());

  std::transform(std::cbegin(L), std::cend(L), std::back_inserter(b_hard),
                 [](const T &e) { return e < 0 ? 1 : 0; });

  return b_hard;
}

template <unsigned iterations, typename U = unsigned, typename R, typename Q,
          typename Func_h, typename Func_v>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum__(const matrix<U> &H, const std::vector<Q> &y, Func_h &&hor,
          Func_v &&vert) {

  matrix<Q> q(H.rows(), H.columns());
  matrix<R> r(H.rows(), H.columns());
  std::vector<R> col_sums(y.size());
  std::vector<R> L(y.size());

  for (unsigned iteration = 0; iteration < iterations; iteration++) {
    vertical__(H, y, r, q, vert);
    horizontal__(H, q, r, hor);

    /* TODO re-use in vertical__ */
    col_sums = column_sum(H, r);

    std::transform(std::cbegin(col_sums), std::cend(col_sums), std::cbegin(y),
                   std::begin(L),
                   [](const R &lhs, const Q &rhs) { return lhs + R(rhs); });
    std::vector<U> b(hard_decision<U>(L));

    if (syndrome(H, b))
      return std::make_tuple(b, L, iteration);
  }

#if 0
  std::cout << H << std::endl;
  
  std::vector<U> b(hard_decision<U>(L));
  for (const auto &bit : b)
    std::cout << static_cast<unsigned>(bit);
  std::cout << std::endl;

  for (const auto &l : L)
    std::cout << l;
  std::cout << std::endl;
#endif
  throw decoding_failure("Decoding failure");
}

template <typename R> R unmodified_horizontal(const R &arg) { return arg; }
template <typename R, typename Q,
          typename Result_t = typename std::common_type<Q, R>::type>
Result_t unmodified_vertical(const R &r, const Q &y, const Q &) {
  return r + R(y);
}

template <typename R> R normalised_horizontal(const R &arg, const R &alpha) {
  return alpha * arg;
}

template <typename R, typename Q>
R normalised_vertical(const R &arg, const Q &y, const Q &, const R &beta) {
  return beta * arg + R(y);
}

template <typename R, typename U = unsigned, typename Q, unsigned Iterations>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        min_sum_tag<Iterations>) {
  return min_sum__<min_sum_tag<Iterations>::iterations, U, R, Q>(
      H, y, unmodified_horizontal<R>, unmodified_vertical<R, Q>);
}

template <typename R, typename U = unsigned, typename Q, unsigned Iterations,
          typename T>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        normalized_min_sum_tag<Iterations, T>) {
  auto alpha = normalized_min_sum_tag<Iterations, T>::alpha;
  return min_sum__<normalized_min_sum_tag<Iterations, T>::iterations, U, R, Q>(
      H, y, std::bind(normalised_horizontal<R>, std::placeholders::_1, alpha),
      unmodified_vertical<R, Q>);
}

template <typename R, typename U = unsigned, typename Q, unsigned Iterations,
          typename T>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        offset_min_sum_tag<Iterations, T>) {
  return min_sum__<offset_min_sum_tag<Iterations, T>::iterations, U, R, Q>(
      H, y, [](const R &min) {
              auto beta = offset_min_sum_tag<Iterations, T>::beta;
              using Result_t =
                  typename std::common_type<R, decltype(beta)>::type;

              return std::max(min - beta, Result_t(0));
            },
      unmodified_vertical<R, Q>);
}

/* http://dud.inf.tu-dresden.de/LDPC/doc/scms/ */
template <typename R, typename U = unsigned, typename Q, unsigned Iterations>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        self_correcting_1_min_sum_tag<Iterations>) {
  return min_sum__<Iterations, U, R, Q>(
      H, y, unmodified_horizontal<R>, [](const R &r, const Q &y_, const Q &q) {
        auto tmp = r + R(y_);
        if (signum(q) == 0 || signum(q) == signum(tmp))
          return tmp;
        else
          return R(0);
      });
}

template <typename R, typename U = unsigned, typename Q, unsigned Iterations>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        self_correcting_2_min_sum_tag<Iterations>) {
  return min_sum__<Iterations, U, R, Q>(
      H, y, unmodified_horizontal<R>, [](const R &r, const Q &y_, const Q &q) {
        auto tmp = r + R(y_);
        if (tmp * R(q) > 0)
          return tmp;
        else
          return R(0.5) * (tmp + R(q));
      });
}

template <typename R, typename U = unsigned, typename Q, unsigned Iterations,
          typename Alpha, typename Beta>
std::tuple<std::vector<U>, std::vector<R>, unsigned>
min_sum(const matrix<U> &H, const std::vector<Q> &y,
        normalized_2d_min_sum_tag<Iterations, Alpha, Beta>) {
  auto alpha = normalized_2d_min_sum_tag<Iterations, Alpha, Beta>::alpha;
  auto beta = normalized_2d_min_sum_tag<Iterations, Alpha, Beta>::beta;
  return min_sum__<Iterations, U, R, Q>(
      H, y, std::bind(normalised_horizontal<R>, std::placeholders::_1, alpha),
      std::bind(normalised_vertical<R, Q>, std::placeholders::_1,
                std::placeholders::_2, std::placeholders::_3, beta));
}
