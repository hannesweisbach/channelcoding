#pragma once

#include <vector>
#include <algorithm>
#include <utility>

namespace gf {

template <typename GF, typename Element = typename GF::element_t>
class polynomial : public std::vector<Element> {
  /* return quotient and remainder */
  std::pair<polynomial, polynomial> division(const polynomial &lhs,
                                             const polynomial &rhs) const {
    if (!rhs)
      throw std::logic_error("Division by zero");

    if (lhs.degree() < rhs.degree())
      return std::make_tuple(polynomial(1, GF::zero), lhs);

    const size_t quotient_degree =
        static_cast<size_t>(lhs.degree() - rhs.degree());
    polynomial q(quotient_degree + 1, GF::zero);
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
                     [=](const Element &lhs_, const Element &rhs_) {
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
                  Element(0));
    }

    std::transform(std::cbegin(rhs), std::cend(rhs), std::cbegin(*this),
                   std::begin(*this), std::forward<Op &&>(op));

    /* todo remove leading zeroes */

    return *this;
  }

  void push_front(const Element &value) {
    this->insert(std::begin(*this), value);
  }

  void push_front(Element &&value) {
    this->insert(std::begin(*this), std::forward<Element &&>(value));
  }

public:
  using std::vector<Element>::vector;
  using element_type = Element;

  //polynomial() : std::vector<Element>(1, GF::zero) {}
  polynomial() = default;
  polynomial(const Element &e) : std::vector<Element>::vector(1, e) {}
  polynomial(const polynomial &) = default;
  polynomial(polynomial &&) = default;
  polynomial &operator=(const polynomial &rhs) = default;
  polynomial &operator=(polynomial &&) = default;

  polynomial(std::vector<typename GF::storage_t> coefficients) {
    for (const auto &element : coefficients) {
      this->push_back(typename GF::element_t(element));
    }
  }
  polynomial(std::vector<int> coefficients) {
    for (const auto &power : coefficients)
      this->push_back(GF::element(power));
  }
  polynomial(std::initializer_list<int> coefficients) {
    for (const auto &power : coefficients)
      this->push_back(GF::element(power));
  }
  explicit polynomial(const std::vector<Element> &v)
      : std::vector<Element>::vector(v) {}

  typename std::vector<Element>::const_reverse_iterator crlead() const
      noexcept {
    return std::find_if(std::crbegin(*this), std::crend(*this),
                        [](const Element &e) { return e; });
  }

  typename std::vector<Element>::reverse_iterator rlead() {
    auto lead = std::find_if(std::rbegin(*this), std::rend(*this),
                             [](const Element &e) { return e; });
    if (lead == std::rend(*this))
      return std::rend(*this) - 1;

    return lead;
  }

  ssize_t degree() const noexcept {
    auto it = std::find_if(std::crbegin(*this), std::crend(*this),
                           [](const Element &e) { return e; });
    return std::distance(it, std::crend(*this)) - 1;
  }

  typename std::vector<Element>::size_type weight() const noexcept {
    return std::count_if(std::cbegin(*this), std::cend(*this),
                         [](const auto &e) { return e; });
  }

  Element highest() const {
    if (!this->size())
      throw std::out_of_range("Polynomial has no terms.");

    return this->at(degree());
  }

  std::vector<Element> zeroes() const {
    GF field;
    std::vector<Element> zeroes;

    for (const auto &element : field) {
      if ((*this)(element) == GF::zero)
        zeroes.push_back(element);
    }
    return zeroes;
  }

  void set(size_t power, const Element &rhs);

  polynomial &reverse() {
    std::reverse(std::begin(*this), std::end(*this));
    return *this;
  }

  const Element &operator[](const size_t &index) const {
    return this->at(index);
  }
  Element &operator[](const size_t &index) { return this->at(index); }

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
      return polynomial({ Element(0) });
    size_t result_degree = static_cast<size_t>(degree() + rhs.degree());
    polynomial result(result_degree + 1, GF::zero);

    auto coffset = std::cbegin(result);
    auto offset = std::begin(result);
    for (const auto &term : *this) {
      auto copy = rhs * term;
      //std::cout << "  rhs * term = " << rhs << " * " << term << " = " << copy << std::endl;
      std::transform(
          std::cbegin(copy), std::cend(copy), coffset, offset,
          [](const Element &lhs, const Element &rhs) { return lhs + rhs; });
      ++coffset;
      ++offset;
      //std::cout << "  " << result << std::endl;
    }
    return result;
  }

  polynomial operator*(const Element &rhs) const {
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
    return element_wise(
        rhs, [](const Element &lhs, const Element &rhs) { return lhs + rhs; });
  }

  polynomial &operator-=(const polynomial &rhs) {
    return element_wise(
        rhs, [](const Element &lhs, const Element &rhs) { return lhs - rhs; });
  }

  polynomial &operator*=(const polynomial &rhs) {
    *this = *this * rhs;
    return *this;
  }

  polynomial &operator*=(const Element &rhs) {
    for (auto &&term : *this)
      term *= rhs;
    return *this;
  }

  friend polynomial operator*(const Element &lhs, const polynomial &rhs) {
    return rhs * lhs;
  }
  
  bool operator==(const polynomial &rhs) {
    return this->operator==(rhs);
  }

  bool operator!=(const polynomial &rhs) {
    return this->operator!=(rhs);
  }

  Element operator()(const Element &x) const {
    if (this->empty())
      return Element(0);

    Element result(*std::crbegin(*this));

    std::for_each(std::crbegin(*this) + 1, std::crend(*this),
                  [&](const Element &coefficient) {
      result = (result * x) + coefficient;
    });
    return result;
  }

  /* returns true for non-empty, non-zero polynomials */
  explicit operator bool() const {
    return this->degree() >= 0;
  }

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
  static polynomial n(unsigned power) {
    polynomial x(power + 1, GF::zero);
    x.back() = GF::one;
    return x;
  }
};

template <typename GF, typename Element>
const polynomial<GF, Element> polynomial<GF, Element>::x = { GF::zero, GF::one };
}
