#pragma once

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <initializer_list>

namespace math {

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

template <typename Coefficient> class polynomial {
  using rep_type = std::vector<Coefficient>;
  using iterator = typename rep_type::iterator;
  using const_iterator = typename rep_type::const_iterator;
  using reverse_iterator = typename rep_type::reverse_iterator;
  using const_reverse_iterator = typename rep_type::const_reverse_iterator;

  rep_type data;

  /* return quotient and remainder */
  std::pair<polynomial, polynomial> division(const polynomial &lhs,
                                             const polynomial &rhs) const {
    if (!rhs)
      throw std::logic_error("Division by zero");

    if (lhs.degree() < rhs.degree())
      return std::make_tuple(polynomial(1, Coefficient(0)), lhs);

    const size_t quotient_degree =
        static_cast<size_t>(lhs.degree() - rhs.degree());
    polynomial q(quotient_degree + 1, Coefficient(0));
    polynomial r(lhs);

    auto r_lead = r.rlead();
    auto rhs_lead = rhs.crlead();
    auto q_it = std::rbegin(q);

    for (; r.degree() >= rhs.degree(); ++r_lead, ++q_it) {
      auto coefficient = *r_lead / *rhs_lead;

      // q = q + t;
      (*q_it) += coefficient;

      // r = r - (rhs * t);
      // leading zero coeffs in rhs possible, so use lead instead of crbegin()
      std::transform(r_lead, r_lead + rhs.degree() + 1, rhs_lead, r_lead,
                     [=](const Coefficient &lhs_, const Coefficient &rhs_) {
        return lhs_ - (rhs_ * coefficient);
      });
    }

    /* TODO remove leading zeroes */

    return std::make_tuple(q, r);
  }

  template <typename Op>
  polynomial &element_wise(const polynomial &rhs, Op &&op) {
    if (data.size() < rhs.data.size()) {
      data.reserve(rhs.data.size());
      std::fill_n(std::back_inserter(data), rhs.data.size() - data.size(),
                  Coefficient(0));
    }

    std::transform(std::cbegin(rhs.data), std::cend(rhs.data),
                   std::cbegin(data), std::begin(data),
                   std::forward<Op &&>(op));

    /* todo remove leading zeroes */

    return *this;
  }

  typename std::vector<Coefficient>::const_reverse_iterator crlead() const
      noexcept {
    return std::find_if(std::crbegin(data), std::crend(data),
                        [](const Coefficient &e) { return e; });
  }

  typename std::vector<Coefficient>::reverse_iterator rlead() {
    auto lead = std::find_if(std::rbegin(data), std::rend(data),
                             [](const Coefficient &e) { return e; });
    if (lead == std::rend(data))
      return std::rend(data) - 1;

    return lead;
  }

public:
  using size_type = typename rep_type::size_type;
  using coefficient_type = Coefficient;
  using value_type = Coefficient;

  polynomial() = default;
  polynomial(const Coefficient &e) : data(1, e) {}
  polynomial(const polynomial &) = default;
  polynomial(std::initializer_list<Coefficient> il) : data(std::move(il)) {}
  polynomial(size_type size, Coefficient c) : data(size, c) {}
  polynomial(polynomial &&) = default;
  polynomial &operator=(const polynomial &rhs) = default;
  polynomial &operator=(polynomial &&) = default;

  explicit polynomial(
      std::vector<typename Coefficient::storage_type> coefficients) {
    for (const auto &element : coefficients) {
      data.push_back(Coefficient(element));
    }
  }
  explicit polynomial(const std::vector<Coefficient> &v) : data(v) {}

  ssize_t degree() const noexcept {
    auto it = std::find_if(std::crbegin(data), std::crend(data),
                           [](const Coefficient &e) { return e; });
    return std::distance(it, std::crend(data)) - 1;
  }

  typename std::vector<Coefficient>::size_type weight() const noexcept {
    return std::count_if(std::cbegin(data), std::cend(data),
                         [](const auto &e) { return e; });
  }

  size_type size() const { return data.size(); }

  iterator begin() noexcept { return data.begin(); }
  iterator end() noexcept { return data.end(); }
  const_iterator begin() const noexcept { return data.begin(); }
  const_iterator end() const noexcept { return data.end(); }
  const_iterator cbegin() const noexcept { return data.cbegin(); }
  const_iterator cend() const noexcept { return data.cend(); }
  reverse_iterator rbegin() noexcept { return data.rbegin(); }
  reverse_iterator rend() noexcept { return data.rend(); }
  const_reverse_iterator rbegin() const noexcept { return data.rbegin(); }
  const_reverse_iterator rend() const noexcept { return data.rend(); }
  const_reverse_iterator crbegin() const noexcept { return data.crbegin(); }
  const_reverse_iterator crend() const noexcept { return data.crend(); }

  void reserve(size_type new_capacity) { data.reserve(new_capacity); }

  void push_back(const value_type &v) { data.push_back(v); }
  void push_back(value_type &&v) { data.push_back(std::move(v)); }

  value_type &at(size_type i) { return data.at(i); }
  const value_type &at(size_type i) const { return data.at(i); }

  rep_type to_vector() const {
    return data;
  }

  Coefficient highest() const {
    if (data.empty())
      throw std::out_of_range("Polynomial has no terms.");

    return data.at(degree());
  }

  polynomial &reverse() {
    std::reverse(std::begin(data), std::begin(data) + degree() + 1);
    return *this;
  }

  polynomial simplified() const {
    polynomial copy;
    copy.data.reserve(data.size());
    std::copy(std::cbegin(data), std::cbegin(data) + this->degree() + 1,
              std::back_inserter(copy.data));

    return copy;
  }

  const Coefficient &operator[](const size_t &index) const {
    return data.at(index);
  }
  Coefficient &operator[](const size_t &index) { return data.at(index); }

  polynomial operator+(const polynomial &rhs) const {
    auto copy(*this);
    copy += rhs;
    return copy;
  }

  polynomial operator-(const polynomial &rhs) const {
    auto copy(*this);
    copy += rhs;
    return copy;
  }

  polynomial operator*(const polynomial &rhs) const {
    if (!(*this && rhs))
      return polynomial({ Coefficient(0) });

    auto rhs_ = rhs.simplified();
    const size_t size = static_cast<size_t>(this->degree()) + rhs_.data.size();
    polynomial result(size, Coefficient(0));

    auto coffset = std::cbegin(result);
    auto offset = std::begin(result);
    for (const auto &term : data) {
      if (term) {
        auto copy = rhs_ * term;

        // std::cout << "  rhs * term = " << rhs << " * " << term << " = " <<
        // copy
        // << std::endl;
        std::transform(std::cbegin(copy), std::cend(copy), coffset, offset,
                       std::plus<Coefficient>{});
      }
      ++coffset;
      ++offset;
    }
    return result;
  }

  polynomial operator*(const Coefficient &rhs) const {
    auto result(*this);
    result *= rhs;
    return result;
  }

  polynomial operator/(const polynomial &rhs) const {
    return division(*this, rhs).first;
  }

  polynomial operator%(const polynomial &rhs) const {
    return division(*this, rhs).second;
  }

  polynomial &operator+=(const polynomial &rhs) {
    return element_wise(rhs, std::plus<Coefficient>{});
  }

  polynomial &operator-=(const polynomial &rhs) {
    return element_wise(rhs, std::minus<Coefficient>{});
  }

  polynomial &operator*=(const polynomial &rhs) {
    *this = *this * rhs;
    return *this;
  }

  polynomial &operator*=(const Coefficient &rhs) {
    for (auto &&term : data)
      term *= rhs;
    return *this;
  }

  friend polynomial operator*(const Coefficient &lhs, const polynomial &rhs) {
    return rhs * lhs;
  }

  bool operator==(const polynomial &rhs) { return data == rhs.data; }
  bool operator!=(const polynomial &rhs) { return data != rhs.data; }

  template <typename T> T operator()(const T &x_) const {
    if (data.empty() || x_ == T(0))
      return T(0);

    T result(*std::crbegin(data));

    std::for_each(std::crbegin(data) + 1, std::crend(data),
                  [&](const Coefficient &coefficient) {
      result = (result * x_) + T(coefficient);
    });
    return result;
  }

  /* returns true for non-empty, non-zero polynomials */
  explicit operator bool() const { return this->degree() >= 0; }

  polynomial gcd(const polynomial &other) const {
    if (!other)
      return *this;
    else
      return other.gcd(*this % other);
  }

  polynomial lcm(const polynomial &other) const {
    if (!*this && !other)
      return polynomial(Coefficient(0));
    else
      return *this / this->gcd(other) * other;
  }

  friend std::ostream &operator<<(std::ostream &s,
                                  const polynomial &polynomial) {

    for (const auto &term : polynomial)
      s << term << " + ";
    return s;
  }

  static const polynomial x;

  /* returns polynomial x^n */
  static polynomial n(size_t power) {
    polynomial p(power + 1, Coefficient(0));
    p.data.back() = Coefficient(1);
    return p;
  }

  template <typename PolyA> static polynomial from(const PolyA &a) {
    polynomial b;
    b.data.reserve(a.data.size());

    std::transform(std::cbegin(a.data), std::cend(a.data),
                   std::back_inserter(b.data),
                   [](const auto &e) { return value_type(e); });

    return b;
  }
};

template <typename Coefficient>
const polynomial<Coefficient> polynomial<Coefficient>::x = { Coefficient(0),
                                                             Coefficient(1) };
}
