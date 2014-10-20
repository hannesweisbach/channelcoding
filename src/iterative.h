#pragma once

#include <vector>
#include <utility>
#include <limits>
#include <numeric>
#include <cmath>
#include <tuple>
#include <functional>

#include "matrix.h"
#include "util.h"

#ifdef NDEBUG
#define at(x) operator[](x)
#endif

template <typename T> constexpr int signum(const T &val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
inline auto syndrome(const matrix<T> &H, const std::vector<T> &b) {
  const auto syndromes = H * b;
  for (const auto &s : syndromes)
    if (s % 2)
      return false;
  return true;
}

/* check node update */
template <typename Q, typename R, typename Functor>
void horizontal__(const matrix<int> &H, const matrix<Q> &q, matrix<R> &r,
                  Functor &&fn) {
  const size_t rows = H.rows();
  const size_t cols = H.columns();
  for (auto row = 0; row < rows; row++) {
    for (auto col = 0; col < cols; col++) {
      if (H.at(row).at(col)) {
        int sign = 1;
        Q min = std::numeric_limits<Q>::max();

        for (size_t i = 0; i < cols; i++) {
          if (i != col && H.at(row).at(i)) {
            sign *= signum(q.at(row).at(i));
            min = std::min(min, std::abs(q.at(row).at(i)));
          }
        }
        r.at(row).at(col) = sign * fn(min);
      }
    }
  }
}

template <typename Q, typename R, typename... Args>
void horizontal(const matrix<int> &H, const matrix<Q> &q, matrix<R> &r,
                Args &&... args) {
  horizontal__<Q, R>(H, q, r, [](const R &arg) { return arg; });
}

template <typename Q, typename R, typename Param>
void horizontal_normalized(const matrix<int> &H, const matrix<Q> &q,
                           matrix<R> &r, const Param &p) {
  horizontal__<Q, R>(H, q, r, [=](const R &min) { return p.alpha() * min; });
}

template <typename Q, typename R, typename Param,
          typename Result_t =
              decltype(std::declval<R>() - std::declval<float>())>
void horizontal_offset(const matrix<int> &H, const matrix<Q> &q, matrix<R> &r,
                       const Param &p) {
  horizontal__(H, q, r, [=](const R &min) {
    return std::max(min - p.beta(), Result_t(0));
  });
}

/* symbol node update */
template <typename Q, typename R, typename Functor>
void vertical__(const matrix<int> &H, const std::vector<Q> y,
                const matrix<R> &r, matrix<Q> &q, Functor &&fn) {
  std::vector<R> col_sums(H.columns(), R(0));
  for (size_t row = 0; row < H.rows(); row++) {
    for (size_t col = 0; col < H.columns(); col++) {
      if (H.at(row).at(col)) {
        col_sums.at(col) += r.at(row).at(col);
      }
    }
  }
  const size_t rows = H.rows();
  const size_t cols = H.columns();
  for (auto row = 0; row < rows; row++) {
    for (auto col = 0; col < cols; col++) {
      if (H.at(row).at(col)) {
        const auto exclusive_colsum = col_sums.at(col) - r.at(row).at(col);
        q.at(row).at(col) = fn(exclusive_colsum, y.at(col), q.at(row).at(col));
      }
    }
  }
}

template <typename Q, typename R, typename... Args>
void vertical(const matrix<int> &H, const std::vector<Q> y, const matrix<R> &r,
              matrix<Q> &q, Args &&...) {
  vertical__<Q, R>(H, y, r, q,
                   [](const R &r, const Q &y, const Q &) { return r + y; });
}

/* http://dud.inf.tu-dresden.de/LDPC/doc/scms/ */
template <typename Q, typename R, typename... Args>
void vertical_sc_1(const matrix<int> &H, const std::vector<Q> y,
                   const matrix<R> &r, matrix<Q> &q, Args &&...) {
  vertical__<Q, R>(H, y, r, q, [](const R &r, const Q &y, const Q &q) {
    auto tmp = r + y;
    if (signum(q) == 0 || signum(q) == signum(tmp))
      return tmp;
    else
      return R(0);
  });
}

template <typename Q, typename R, typename Param>
void vertical_normalized(const matrix<int> &H, const std::vector<Q> y,
                         const matrix<R> &r, matrix<Q> &q, const Param &p) {
  vertical__<Q, R>(H, y, r, q, [=](const R &r, const Q &y,
                                   const Q &) { return r * p.beta() + y; });
}

template <typename Q, typename R>
void vertical_sc_2(const matrix<int> &H, const std::vector<Q> y,
                   const matrix<R> &r, matrix<Q> &q) {
  vertical__<Q, R>(H, y, r, q, [](const R &r, const Q &y, const Q &q) {
    auto tmp = r + y;
    if (tmp * q > 0)
      return tmp;
    else
      return R(0.5) * (tmp + q);
  });
}

template <typename Q, typename R,
          typename Result_t = decltype(std::declval<Q>() * std::declval<R>())>
std::vector<Result_t> likelihood(const std::vector<Q> &y, const matrix<int> &H,
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

template <typename T> std::vector<int> hard_decision(const std::vector<T> &L) {
  std::vector<int> b_hard;
  b_hard.reserve(L.size());

  std::transform(std::cbegin(L), std::cend(L), std::back_inserter(b_hard),
                 [](const T &e) { return e < 0 ? 1 : 0; });

  return b_hard;
}

template <unsigned iterations, typename R, typename Q, typename Func_h,
          typename Func_v, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
min_sum__(const matrix<int> &H, const std::vector<Q> &y, Func_h hor,
          Func_v vert, Args &&... args) {

  matrix<Q> q(H.rows(), H.columns());
  matrix<R> r(H.rows(), H.columns());
  std::vector<R> L;
  std::copy(std::cbegin(y), std::cend(y), std::back_inserter(L));
  std::vector<std::vector<int>> history;

  for (unsigned iteration = 0; iteration < iterations; iteration++) {
    vert(H, y, r, q, std::forward<Args>(args)...);
    hor(H, q, r, std::forward<Args>(args)...);
#if 0
    vertical(H, L, r, q);
    horizontal(H, q, r);
    //horizontal_normalized(H, q, r, 1.0f);
    //horizontal_offset(H, q, r, 1.0f);
#endif
#if 0
    std::cout << "Q: " << std::endl << q << std::endl;
    std::cout << "R: " << std::endl << r << std::endl;
#endif
    
    L = likelihood(y, H, r);
    std::vector<int> b(hard_decision(L));
    
#if 0
    for (auto e : L)
      std::cout << std::setw(2) << e << " ";
    std::cout << std::endl;
#endif
    
#if 0
    for (auto e : b)
      std::cout << std::setw(2) << e << " ";
    std::cout << std::endl;
#endif
    //history.push_back(b);

    if (syndrome(H, b))
      return std::make_tuple(b, L, iteration);
  }
#if 0
  for (const auto &b_h : history) {
    for (auto e : b_h)
      std::cout << std::setw(2) << e << " ";
    std::cout << std::endl;
  }
#endif
  throw std::runtime_error("Decoding failure");
}

template <unsigned iterations, typename R, typename Q, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
min_sum(const matrix<int> &H, const std::vector<Q> &y) {
  return min_sum__<iterations, R, Q>(H, y, horizontal<Q, R, Args...>,
                                     vertical<Q, R>);
}

class alpha_wrapper {
  const float alpha_;

public:
  alpha_wrapper(const float alpha) : alpha_(alpha) {}
  float alpha() const { return alpha_; }
};

class beta_wrapper {
  const float beta_;

public:
  beta_wrapper(const float beta) : beta_(beta) {}
  float beta() const { return beta_; }
};

class ab_wrapper {
  const float alpha_;
  const float beta_;

public:
  ab_wrapper(const float alpha, const float beta)
      : alpha_(alpha), beta_(beta) {}
  float alpha() const { return alpha_; }
  float beta() const { return beta_; }
};

template <unsigned iterations, typename R, typename Q, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
nms(const matrix<int> &H, const std::vector<Q> &y, const float alpha = 1.0f) {
  return min_sum__<iterations, R, Q>(
      H, y, horizontal_normalized<Q, R, alpha_wrapper>,
      vertical<Q, R, alpha_wrapper>, alpha_wrapper(alpha));
}

template <unsigned iterations, typename R, typename Q, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
oms(const matrix<int> &H, const std::vector<Q> &y, const float beta = 0.0f) {
  return min_sum__<iterations, R, Q>(
      H, y, horizontal_offset<Q, R, beta_wrapper>, vertical<Q, R, beta_wrapper>,
      beta_wrapper(beta));
}

template <unsigned iterations, typename R, typename Q, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
scms1(const matrix<int> &H, const std::vector<Q> &y) {
  return min_sum__<iterations, R, Q>(H, y, horizontal<Q, R>,
                                     vertical_sc_1<Q, R>);
}

template <unsigned iterations, typename R, typename Q, typename... Args>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
scms2(const matrix<int> &H, const std::vector<Q> &y) {
  return min_sum__<iterations, R, Q>(H, y, horizontal<Q, R>,
                                     vertical_sc_2<Q, R>);
}

template <unsigned iterations, typename R, typename Q>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
nms_2d(const matrix<int> &H, const std::vector<Q> &y, const float alpha = 1.0f,
       const float beta = 0.1f) {
  return min_sum__<iterations, R, Q>(
      H, y, horizontal_normalized<Q, R, ab_wrapper>,
      vertical_normalized<Q, R, ab_wrapper>, ab_wrapper(alpha, beta));
}

template <unsigned iterations, typename R, typename Q>
std::tuple<std::vector<int>, std::vector<R>, unsigned>
scms1_nms(const matrix<int> &H, const std::vector<Q> &y,
          const float alpha = 1.0f) {
  return min_sum__<iterations, R, Q>(
      H, y, horizontal_normalized<Q, R, alpha_wrapper>,
      vertical_sc_1<Q, R, alpha_wrapper>, alpha_wrapper(alpha));
}
