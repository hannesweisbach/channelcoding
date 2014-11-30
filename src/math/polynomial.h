#pragma once

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>

namespace math {

template <typename Coefficient>
class polynomial : public std::vector<Coefficient> {
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
    if (this->size() < rhs.size()) {
      this->reserve(rhs.size());
      std::fill_n(std::back_inserter(*this), rhs.size() - this->size(),
                  Coefficient(0));
    }

    std::transform(std::cbegin(rhs), std::cend(rhs), std::cbegin(*this),
                   std::begin(*this), std::forward<Op &&>(op));

    /* todo remove leading zeroes */

    return *this;
  }

  typename std::vector<Coefficient>::const_reverse_iterator crlead() const
      noexcept {
    return std::find_if(std::crbegin(*this), std::crend(*this),
                        [](const Coefficient &e) { return e; });
  }

  typename std::vector<Coefficient>::reverse_iterator rlead() {
    auto lead = std::find_if(std::rbegin(*this), std::rend(*this),
                             [](const Coefficient &e) { return e; });
    if (lead == std::rend(*this))
      return std::rend(*this) - 1;

    return lead;
  }

public:
  using std::vector<Coefficient>::vector;
  using coefficient_type = Coefficient;

  polynomial() = default;
  polynomial(const Coefficient &e) : std::vector<Coefficient>::vector(1, e) {}
  polynomial(const polynomial &) = default;
  polynomial(polynomial &&) = default;
  polynomial &operator=(const polynomial &rhs) = default;
  polynomial &operator=(polynomial &&) = default;

  explicit polynomial(
      std::vector<typename Coefficient::storage_t> coefficients) {
    for (const auto &element : coefficients) {
      this->push_back(Coefficient(element));
    }
  }
  explicit polynomial(const std::vector<Coefficient> &v)
      : std::vector<Coefficient>::vector(v) {}

  ssize_t degree() const noexcept {
    auto it = std::find_if(std::crbegin(*this), std::crend(*this),
                           [](const Coefficient &e) { return e; });
    return std::distance(it, std::crend(*this)) - 1;
  }

  typename std::vector<Coefficient>::size_type weight() const noexcept {
    return std::count_if(std::cbegin(*this), std::cend(*this),
                         [](const auto &e) { return e; });
  }

  Coefficient highest() const {
    if (!this->size())
      throw std::out_of_range("Polynomial has no terms.");

    return this->at(degree());
  }

  polynomial &reverse() {
    std::reverse(std::begin(*this), std::begin(*this) + degree() + 1);
    return *this;
  }

  const Coefficient &operator[](const size_t &index) const {
    return this->at(index);
  }
  Coefficient &operator[](const size_t &index) { return this->at(index); }

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
    /* The range based for loop below also loops over leading zero coefficients,
     * so use size() instead of degree()
     */
    polynomial result(this->size() + rhs.size(), Coefficient(0));

    auto coffset = std::cbegin(result);
    auto offset = std::begin(result);
    for (const auto &term : *this) {
      auto copy = rhs * term;

      // std::cout << "  rhs * term = " << rhs << " * " << term << " = " << copy
      // << std::endl;
      std::transform(std::cbegin(copy), std::cend(copy), coffset, offset,
                     std::plus<Coefficient>{});
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
    for (auto &&term : *this)
      term *= rhs;
    return *this;
  }

  friend polynomial operator*(const Coefficient &lhs, const polynomial &rhs) {
    return rhs * lhs;
  }

  bool operator==(const polynomial &rhs) { return this->operator==(rhs); }

  bool operator!=(const polynomial &rhs) { return this->operator!=(rhs); }

  Coefficient operator()(const Coefficient &x_) const {
    if (this->empty() || x_ == Coefficient(0))
      return Coefficient(0);

    Coefficient result(*std::crbegin(*this));

    std::for_each(std::crbegin(*this) + 1, std::crend(*this),
                  [&](const Coefficient &coefficient) {
      result = (result * x_) + coefficient;
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
      return polynomial(0);
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
    p.back() = Coefficient(1);
    return p;
  }
};

template <typename Coefficient>
const polynomial<Coefficient> polynomial<Coefficient>::x = { Coefficient(0),
                                                             Coefficient(1) };
}
