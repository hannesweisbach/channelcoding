#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <numeric>

#include "util.h"
#include "gf.h"

gf::gf(unsigned q, uint64_t modular_polynomial)
    : log([=]() {
        const unsigned n = (1 << q) - 1;

        if (!(modular_polynomial & (1 << q))) {
          std::ostringstream os;
          os << "The " << q << "-th term is not 1 in the modular polynomial "
             << std::hex << modular_polynomial << std::dec << std::endl;
          throw std::runtime_error(os.str());
        }
        if ((modular_polynomial & n) > (1 << q)) {
          std::ostringstream os;
          os << "The " << q << "-th term is not the highest non-zero "
                               "term in moduluar polynomial " << std::hex
             << modular_polynomial << std::dec << std::endl;

          throw std::runtime_error(os.str());
        }
        std::vector<uint8_t> table(2 * (n + 1), -1);
        uint8_t power = 0;
        uint8_t polynomial = 1;
        for (unsigned power = 0; power < n; power++) {
          table.at(polynomial) = power;
          table.at(polynomial + n + 1) = power;
          
          bool carry = polynomial & (1 << (q - 1));
          polynomial <<= 1;
          if (carry)
            polynomial ^= modular_polynomial;
        }
        return table;
      }()),
      exp([=]() {
        const unsigned n = (1 << q) - 1;
        std::vector<uint8_t> table(2 * n, 0);
        uint8_t power = 0;
        uint8_t polynomial = 1;
        for (unsigned power = 0; power < n; power++) {
          table.at(power) = polynomial;
          table.at(power + n) = polynomial;

          bool carry = polynomial & (1 << (q - 1));
          polynomial <<= 1;
          if (carry)
            polynomial ^= modular_polynomial;
        }
        return table;
      }()) {}

gf::const_iterator gf::begin() const { return gf_element(*this, 0); }
gf::const_iterator gf::end() const {
  return gf_element(*this, exp.size() / 2 + 1);
}

size_t gf::size() const { return exp.size() / 2; }

void gf::dump() const {
  for (unsigned alpha = 0; alpha < log.size() / 2; alpha++)
    std::cout << "alpha: " << std::hex << std::setw(2) << alpha << std::dec
              << " " << std::setw(2) << (int)log.at(alpha) << std::endl;
  
  std::cout << std::endl;

  for (unsigned power = 0; power < exp.size() / 1; power++)
    std::cout << "α^" << power << " = "
              << std::hex << static_cast<unsigned>(exp.at(power)) << std::dec
              << std::endl;
}

gf_element gf::power_to_polynomial(const unsigned &power) const {
  return gf_element(*this, exp.at(power));
}

unsigned gf::polynomial_to_power(const gf_element &element) const {
  return log.at(element.value);
}

gf_element gf::zero() const { return gf_element(*this, 0); }
gf_element gf::one() const { return gf_element(*this, 1); }
gf_element gf::alpha() const { return gf_element(*this, 2); }

gf::const_iterator::const_iterator(const gf_element &e) : e(e) {}

const gf_element &gf::const_iterator::operator*() const { return e; }
const gf_element *gf::const_iterator::operator->() const { return &e; }

const gf::const_iterator &gf::const_iterator::operator++() {
  ++e.value;
  return *this;
}

gf::const_iterator gf::const_iterator::operator++(int) {
  const_iterator result(*this);
  ++(*this);
  return result;
}

bool gf::const_iterator::operator!=(const const_iterator &rhs) {
  return e != *rhs;
}

bool gf::operator==(const gf &rhs) const {
  return std::equal(std::begin(log), std::end(log), std::begin(rhs.log),
                    std::end(rhs.log)) &&
         std::equal(std::begin(exp), std::end(exp), std::begin(rhs.exp),
                    std::end(rhs.exp));
}

bool gf::operator!=(const gf &rhs) const { return !(*this == rhs); }

#if 0
gf_element gf::pow(const gf_element &lhs, const unsigned power) const {
  return gf_element(*this, exp.at((log.at(lhs.value) * power) % exp.size()));
}
#endif

gf_element::gf_element(const gf &field, uint8_t value)
    : field(&field), value(value) {}

gf_element gf_element::operator+(const gf_element &rhs) const {
  assert(rhs.field);
  return gf_element(*field, value ^ rhs.value);
}

void gf_element::assert(const gf *const other) const {
  if(*field != *other) {
    std::cout << "This field: " << field << std::endl;
    field->dump();
    std::cout << std::endl << "Other field: " << other << std::endl;
    other->dump();
    throw std::runtime_error("Incompatible gf elements");
  }
}


gf_element gf_element::operator-(const gf_element &rhs) const {
  assert(rhs.field);
  return *this + rhs;
}

gf_element gf_element::operator*(const gf_element &rhs) const {
  assert(rhs.field);
  if (value == 0 || rhs.value == 0)
    return gf_element(*field, 0);
  return field->power_to_polynomial(field->polynomial_to_power(*this) +
                                    field->polynomial_to_power(rhs));
}

gf_element gf_element::operator/(const gf_element &rhs) const {
  assert(rhs.field);
  if(rhs.value == 0)
    throw std::overflow_error("Divide by zero exception");

  if (value == 0)
    return gf_element(*field, 0);

  return field->power_to_polynomial(field->polynomial_to_power(*this) -
                                    field->polynomial_to_power(rhs) +
                                    field->size());
}

gf_element &gf_element::operator*=(const gf_element &rhs) {
  *this = *this * rhs;
  return *this;
}
  
gf_element &gf_element::operator++() {
  if(value == 0) {
    value = 1;
  } else {
    *this = field->power_to_polynomial(power() + 1);
  }
  return *this;
}

gf_element gf_element::operator++(int) {
  gf_element result(*this);
  ++(*this);
  return result;
}

bool gf_element::operator!=(const gf_element &rhs) const {
  assert(rhs.field);
  return value != rhs.value;
}

bool gf_element::operator==(const gf_element &rhs) const {
  assert(rhs.field);
  return value == rhs.value;
}

gf_element gf_element::inverse() const {
  return field->power_to_polynomial(0) / *this;
}

gf_element gf_element::zero() const { return gf_element(*field, 0); }

gf_element::operator bool() const { return value != 0; }

unsigned gf_element::power() const { return field->polynomial_to_power(*this); }

std::ostream &operator<<(std::ostream &os, const gf_element &e) {
  if(e.value == 0)
    return os << 0;
  else
    return os << "α^" << e.power();
}

gf_polynomial::gf_polynomial(const gf &field) : field(&field) {}

gf_polynomial::gf_polynomial(const gf &field,
                             std::vector<gf_element> coefficients)
    : std::vector<gf_element>(std::move(coefficients)), field(&field) {}

gf_polynomial::gf_polynomial(const gf &field,
                             std::initializer_list<gf_element> list)
    : std::vector<gf_element>(std::move(list)), field(&field) {}

gf_polynomial::gf_polynomial(const gf &field, std::vector<int> coefficients,
                             bool powers)
    : field(&field) {
  if (powers) {
    for (const auto &power : coefficients)
      push_back(field.power_to_polynomial(power));
  } else {
    for (const auto &element : coefficients)
      emplace_back(field, element);
  }
}

unsigned gf_polynomial::degree() const {
  auto it = std::find_if(std::crbegin(*this), std::crend(*this),
                         [](const gf_element &e) { return e; });
  auto distance = std::distance(it, std::crend(*this));
  if (distance > 0)
    distance--;

  return distance;
}

const gf_element &gf_polynomial::operator[](const size_t &index) const {
  return at(index);
}
gf_element &gf_polynomial::operator[](const size_t &index) { return at(index); }

void gf_polynomial::push_front(const gf_element &value) {
  insert(std::begin(*this), value);
}

void gf_polynomial::push_front(gf_element &&value) {
  insert(std::begin(*this), std::move(value));
}

gf_polynomial &gf_polynomial::reverse() {
  std::reverse(std::begin(*this), std::end(*this));
  return *this;
}

void gf_polynomial::set(size_t power, const gf_element &rhs) {
  if (power < size()) {
    at(power) = rhs;
  } else {
    insert(std::end(*this), power - size(), rhs.zero());
    push_back(rhs);
  }
}

unsigned gf_polynomial::weight() const {
  return std::count_if(std::cbegin(*this), std::cend(*this),
                       [](const auto &e) { return e; });
}

gf_element gf_polynomial::highest() const {
  auto index = degree();
  if (index) {
    return at(index);
  } else if (size()) {
    return at(0);
  } else {
    return field->zero();
  }
}

std::vector<gf_element> gf_polynomial::zeroes() const {
  std::vector<gf_element> zeroes;

  for (const auto &element : *field) {
    if ((*this)(element) == field->zero())
      zeroes.push_back(element);
  }
  return zeroes;
}

gf_polynomial &gf_polynomial::operator=(const gf_polynomial &rhs) = default;

gf_polynomial gf_polynomial::operator+(const gf_polynomial &rhs) const {
  auto minmax = std::minmax(*this, rhs, [](const auto &lhs, const auto &rhs) {
    return lhs.size() < rhs.size();
  });
  const auto &smaller = minmax.first;
  auto tmp(minmax.second);

  std::transform(
      smaller.cbegin(), smaller.cend(), tmp.cbegin(), std::begin(tmp),
      [](const gf_element &lhs, const gf_element &rhs) { return rhs + lhs; });

  return tmp;
}

gf_polynomial &gf_polynomial::operator+=(const gf_polynomial &rhs) {
  *this = *this + rhs;
  return *this;
}

gf_polynomial gf_polynomial::operator-(const gf_polynomial &rhs) const {
  return *this + rhs;
}

gf_polynomial &gf_polynomial::operator-=(const gf_polynomial &rhs) {
  *this = *this - rhs;
  return *this;
}

gf_polynomial gf_polynomial::operator*(const gf_polynomial &rhs) const {
  gf_polynomial copy(rhs);
  gf_polynomial result(*field, { field->zero() });
  
  for (const auto &term : *this) {
    result += copy * term;
    copy.push_front(field->zero());
  }
  return result;
}

gf_polynomial &gf_polynomial::operator*=(const gf_polynomial &rhs) {
  *this = *this * rhs;
  return *this;
}

gf_polynomial gf_polynomial::operator*(const gf_element &rhs) const {
  auto tmp = *this;
  for (auto &&coefficient : tmp)
    coefficient *= rhs;
  return tmp;
}

gf_polynomial &gf_polynomial::operator*=(const gf_element &rhs) {
  *this = *this * rhs;
  return *this;
}

void gf_polynomial::division(const gf_polynomial &lhs, const gf_polynomial &rhs,
                             gf_polynomial &q, gf_polynomial &r) const {
  r = lhs;
  q.clear();

  for (; r.degree() && r.degree() >= rhs.degree();) {
    //std::cout << r.degree() << " " << rhs.degree() << std::endl;
    auto coefficient = r.highest() / rhs.highest();
    gf_polynomial t(*field);
    std::fill_n(std::back_inserter(t), r.degree() - rhs.degree() + 1,
                field->zero());
    t.back() = r.highest() / rhs.highest();
    //std::cout << "t = " << t << std::endl;
    //q = q * x +  t;
    q = q + t;
    //std::cout << "q = " << q << std::endl;
    //std::cout << "rhs * q = " << rhs * t << std::endl;
    r = r - (rhs * t);
    //std::cout << "r = " << r << std::endl;
  }
  //std::cout << " = " << q << " (rest " << r << ")" << std::endl;
}

gf_polynomial gf_polynomial::operator/(const gf_polynomial &rhs) const {
  gf_polynomial q(*field);
  gf_polynomial r(*this);
  division(*this, rhs, q, r);
  return q;
}

gf_polynomial gf_polynomial::operator%(const gf_polynomial &rhs) const {
  gf_polynomial q(*field);
  gf_polynomial r(*this);
  division(*this, rhs, q, r);
  return r;
}

bool gf_polynomial::operator==(const gf_polynomial &rhs) {
  return (field == rhs.field) &&
         std::equal(std::cbegin(*this), std::cend(*this), std::cbegin(rhs),
                    std::cend(rhs));
}

bool gf_polynomial::operator!=(const gf_polynomial &rhs) {
  return !(*this == rhs);
}

gf_element gf_polynomial::operator()(const gf_element &x) const {
  gf_element result(*std::crbegin(*this));

  std::for_each(std::crbegin(*this) + 1, std::crend(*this),
                [&](const gf_element &coefficient) {
#if 0
    std::cout << result << " * " << x << " = " << field.mul(result, x) << " + "
              << coefficient << " = " << field.mul(result, x) + coefficient << std::endl;
#endif
    result = (result * x) + coefficient;
  });
  return result;
}

std::ostream &operator<<(std::ostream &s, const gf_polynomial &polynomial) {
  for (const auto &term : polynomial)
    s << term << " + ";
  return s;
}

gf_matrix::gf_matrix(const gf &field) : field(&field) {}
gf_matrix::gf_matrix(const gf &field, std::vector<gf_polynomial> rows)
    : std::vector<gf_polynomial>(std::move(rows)), field(&field) {}

gf_polynomial gf_matrix::get_column(const size_t column) const {
  gf_polynomial col(*field);

  for (const auto &row : *this) {
    col.push_back(row.at(column));
  }
  return col;
}

bool gf_matrix::has_zero_row() const {
  auto zero_row =
      std::find_if(std::cbegin(*this), std::cend(*this), [](const auto &row) {
        for (const auto &element : row)
          if (element)
            return false;
        return true;
      });
  return zero_row != std::cend(*this);
}

bool gf_matrix::has_zero_column() const {
  const size_t rows = size();
  if (!rows)
    return true;
  const size_t columns = at(0).size();
  for (size_t column = 0; column < columns; column++) {
    bool zero = true;
    for (size_t row = 0; row < rows; row++) {
      if (at(row).at(column)) {
        zero = false;
        break;
      }
    }
    if (zero)
      return true;
  }

  return false;
}

gf_matrix gf_matrix::reduced_echelon_form() const {
  std::vector<gf_polynomial> nrows(*this);

  /* make sure rows are sorted with the left-most elements at the top */
  std::sort(std::begin(nrows), std::end(nrows),
            [](const gf_polynomial &lhs, const gf_polynomial &rhs) {
    return lhs.degree() > rhs.degree();
  });

  for (auto first_row = nrows.begin(); first_row != nrows.end(); ++first_row) {
    (*first_row) *= (*first_row)[first_row->degree()].inverse();
    for (auto next_rows = first_row + 1; next_rows != nrows.end();
         ++next_rows) {
      /* only add if not already zero */
      if (next_rows->at(first_row->degree()))
        (*next_rows) += (*first_row) * (*next_rows)[next_rows->degree()];
    }
    // std::cout << gf_matrix(*field, nrows) << std::endl;
  }

  for (auto modify = nrows.rbegin() + 1; modify != nrows.rend(); ++modify) {
    for (auto row = nrows.crbegin(); row != modify; ++row) {
      const size_t index = std::distance(std::crbegin(nrows), row) + 1;
      auto factor = (*modify)[index];
      (*modify) += (*row) * factor;
    }
  }

  return gf_matrix(*field, nrows);
}

gf_polynomial gf_matrix::solution() const {
  gf_matrix copy = reduced_echelon_form();
  
  if (copy.back()[1] == field->zero())
    throw std::runtime_error("Linear equation system not solvable");

  gf_polynomial s(*field);
  for (const auto &row : copy) {
    s.push_back(row[0]);
  }

  return s;
}

size_t gf_matrix::columns() const {
  auto max = std::max_element(
      std::begin(*this), std::end(*this),
      [](const auto &lhs, const auto &rhs) { return lhs.size() < rhs.size(); });
  if (max == std::end(*this))
    return 0;
  else
    return max->size();
}

size_t gf_matrix::rows() const { return size(); }

gf_polynomial gf_matrix::operator*(const gf_polynomial &rhs) const {
  gf_polynomial result(*field);

  for (const auto &row : *this) {
    auto first1 = (row.size() < rhs.size()) ? std::begin(row) : std::begin(rhs);
    auto last1 = (row.size() < rhs.size()) ? std::end(row) : std::end(rhs);
    auto first2 = (row.size() < rhs.size()) ? std::begin(rhs) : std::begin(row);

    result.push_back(std::inner_product(first1, last1, first2, field->zero()));
  }
  return result;
}

std::ostream &operator<<(std::ostream &s, const gf_matrix &matrix) {
  for (const auto &row : matrix)
    s << row << "[" << row.degree() << "]" << std::endl;
  return s;
}

