#pragma once

#include <vector>
#include <cstdint>
#include <iostream>
#include <initializer_list>

class gf;
class gf_polynomial;
class gf_matrix;

class gf_element {
  const gf * field;
  uint8_t value;

  void assert(const gf *const other) const;
  friend class gf;
  friend class gf_polynomial;
  friend class gf_matrix;

public:
  gf_element(const gf &field, uint8_t value);
  
  unsigned power() const;
  
  gf_element operator+(const gf_element &rhs) const;
  gf_element operator-(const gf_element &rhs) const;
  gf_element operator*(const gf_element &rhs) const;
  gf_element operator/(const gf_element &rhs) const;
  gf_element &operator*=(const gf_element &rhs);
  gf_element &operator++();
  gf_element operator++(int);
  bool operator!=(const gf_element &rhs) const;
  bool operator==(const gf_element &rhs) const;
  gf_element inverse() const;
  gf_element zero() const;

  operator bool() const;
  friend std::ostream &operator<<(std::ostream &, const gf_element &);
};


class gf {
  const std::vector<uint8_t> log;
  const std::vector<uint8_t> exp;

public:
  class const_iterator {
    //int power;
    gf_element e;

  public:
    const_iterator(const gf_element &);
    const gf_element &operator*() const;
    const gf_element *operator->() const;
    const const_iterator &operator++();
    const_iterator operator++(int);
    bool operator!=(const const_iterator &rhs);
  };
  gf(unsigned q, uint64_t modular_polynomial);

  const_iterator begin() const;
  const_iterator end() const;

  size_t size() const;
  gf_element power_to_polynomial(const unsigned &element) const;
  unsigned polynomial_to_power(const gf_element &element) const;
  gf_element zero() const;
  gf_element one() const;
  gf_element alpha() const;
  void dump() const;

  bool operator==(const gf &rhs) const;
  bool operator!=(const gf &rhs) const;
};

class gf_polynomial : public std::vector<gf_element> {
  const gf * field;

  void division(const gf_polynomial &lhs, const gf_polynomial &rhs,
                gf_polynomial &q, gf_polynomial &r) const;

public:
  gf_polynomial(const gf &field);
  gf_polynomial(const gf &field, std::vector<gf_element> coefficients);
  gf_polynomial(const gf &field, std::initializer_list<gf_element> list);
  gf_polynomial(const gf &field, std::vector<int> coefficients,
                bool powers = false);

  unsigned degree() const;
  unsigned weight() const;
  gf_element highest() const;
  std::vector<gf_element> zeroes() const;

  void set(size_t power, const gf_element &rhs);
  void push_front(const gf_element &value);
  void push_front(gf_element &&value);
  gf_polynomial &reverse();

  const gf_element &operator[](const size_t &index) const;
  gf_element &operator[](const size_t &index);

  gf_polynomial &operator=(const gf_polynomial &rhs);
  gf_polynomial operator+(const gf_polynomial &rhs) const;
  gf_polynomial &operator+=(const gf_polynomial &rhs);
  gf_polynomial operator-(const gf_polynomial &rhs) const;
  gf_polynomial &operator-=(const gf_polynomial &rhs);
  gf_polynomial operator*(const gf_polynomial &rhs) const;
  gf_polynomial &operator*=(const gf_polynomial &rhs);
  gf_polynomial operator*(const gf_element &rhs) const;
  gf_polynomial &operator*=(const gf_element &rhs);
  gf_polynomial operator/(const gf_polynomial &rhs) const;
  gf_polynomial operator%(const gf_polynomial &rhs) const;
  bool operator==(const gf_polynomial &rhs);
  bool operator!=(const gf_polynomial &rhs);
  gf_element operator()(const gf_element &x) const;
  // add();
  // mul();
  friend std::ostream &operator<<(std::ostream &, const gf_polynomial &);
};

//this should rather be gf_linear_equation_system …
class gf_matrix : public std::vector<gf_polynomial> {
  const gf *field;

  gf_polynomial get_column(const size_t column) const;
  bool has_zero_column() const;
  bool has_zero_row() const;
public:
  gf_matrix(const gf &field);
  gf_matrix(const gf &field, std::vector<gf_polynomial> rows);

  gf_matrix reduced_echelon_form() const;
  gf_polynomial solution() const;

  size_t columns() const;
  size_t rows() const;

  gf_polynomial operator*(const gf_polynomial &rhs) const;
  friend std::ostream &operator<<(std::ostream &, const gf_matrix &);
};
