#pragma once

#include <array>

namespace math {
namespace detail {

template <long i> struct log2 {
  static_assert(i > 0, "log2 of 0 cannot be taken.");
  static constexpr long value = log2<i / 2>::value + 1;
};

template <> struct log2<1> {
  static constexpr long value = 0;
};

static constexpr std::array<unsigned, 9> modular_polynomials = {
  { 0, 0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d }
};

/* uint16_t allows galois fields up to 2^15 to be used. */
template <uint16_t poly> struct modular_polynomial {
  static constexpr uint16_t value = poly;
};

template <int size_type> struct uint_t_helper;
template <> struct uint_t_helper<0> {
  using type = uint8_t;
};
template <> struct uint_t_helper<1> {
  using type = uint16_t;
};
template <> struct uint_t_helper<2> {
  using type = uint32_t;
};
template <> struct uint_t_helper<3> {
  using type = uint64_t;
};
template <> struct uint_t_helper<4> {
  using type = __uint128_t;
};

template <unsigned bits> struct uint_t {
  static_assert(bits <= std::numeric_limits<uint64_t>::digits,
                "More than 64 Bits are not supported.");

  using type = typename uint_t_helper<
      (bits > std::numeric_limits<uint8_t>::digits) +
      (bits > std::numeric_limits<uint16_t>::digits) +
      (bits > std::numeric_limits<uint32_t>::digits) +
      (bits > std::numeric_limits<uint64_t>::digits)>::type;
};
}

/* uint16_t allows galois fields up to 2^15 to be used. */
template <uint16_t poly> struct modular_polynomial {
  static constexpr uint16_t value = poly;
};

namespace detail {
template <unsigned q> struct default_modular_polynomial {
  static_assert(q > 0, "GF(2^0) does not make sense. Choose q > 0.");
  static_assert(q < 9, "modular polynomial for GF(q), q > 2^8 have to "
                       "specified manually (or the list of default modular "
                       "polynomials has to be extended.");
  using type =
      ::math::modular_polynomial<std::get<q>(detail::modular_polynomials)>;
};
}

template <long prime, long power,
          typename Modular_Polynomial =
              typename detail::default_modular_polynomial<power>::type>
class ef_element {
  static constexpr size_t size = (1 << power);
  static constexpr unsigned mod_polynomial = Modular_Polynomial::value;
  /* number of non-zero elements in the field */
  static constexpr size_t n = size - 1;
  using storage_type = typename detail::uint_t<power>::type;
};

template <typename GF>
static constexpr std::pair<typename GF::Exp_table_type,
                           typename GF::Log_table_type>
init_tables(const long dimension, const unsigned long size,
            const uint16_t modular_polynomial);

template <long Power, typename Modular_Polynomial>
class ef_element<2, Power, Modular_Polynomial> {
  static constexpr size_t size = (1 << Power);
  /* number of non-zero elements in the field */
  static constexpr size_t n = size - 1;

  static_assert(
      std::is_same<Modular_Polynomial,
                   modular_polynomial<Modular_Polynomial::value> >::value,
      "Modular_Polynomial has to be of type modular_polynomial<>.");
  static_assert(Power <= 64, "Only power <= 64 are supported.");
  static_assert(Power, "power must be non-zero.");
  static_assert(Modular_Polynomial::value & (1 << Power),
                "Power-th bit of modular polynomial has to be set.");
  static_assert(
      (Modular_Polynomial::value & ~(size - 1)) <= size,
      "Power-th bit has to be the highest set bit in the modular polynomial.");

public:
  using storage_type = typename detail::uint_t<Power>::type;

private:
  using element_type = ef_element;
  using Log_table_type = std::array<storage_type, 2 * size>;
  using Exp_table_type = std::array<ef_element, 2 * size>;

  static const Log_table_type log;
  static const Exp_table_type exp;

  template <typename GF>
  friend constexpr std::pair<typename GF::Exp_table_type,
                             typename GF::Log_table_type>
  init_tables(const long dimension, const unsigned long size,
              const uint16_t modular_polynomial);

  struct ef {
    using const_iterator = typename Exp_table_type::const_iterator;
    using const_reverse_iterator =
        typename Exp_table_type::const_reverse_iterator;
    static constexpr size_t offset = size;
    using element_type = ef_element;

    const_iterator begin() const { return std::begin(exp) + offset; }
    const_iterator end() const { return std::end(exp); }
    const_iterator cbegin() const noexcept { return std::cbegin(exp) + offset; }
    const_iterator cend() const noexcept { return std::cend(exp); }
    const_reverse_iterator rbegin() const { return std::rbegin(exp); }
    const_reverse_iterator rend() const { return std::rend(exp) - offset; }
    const_reverse_iterator crbegin() const noexcept {
      return std::crbegin(exp);
    }
    const_reverse_iterator crend() const noexcept {
      return std::crend(exp) - offset;
    }
    friend std::ostream &operator<<(std::ostream &os, const ef &ef) {
      for (const auto &e : exp) {
        os << e << "(" << static_cast<unsigned>(e.value) << ")" << std::endl;
      }
      for (const auto &e : ef) {
        os << e << "(" << static_cast<unsigned>(e.value) << ")" << std::endl;
      }
      return os;
    }
  };

  friend ef;
  storage_type value = 0;

public:
  using field_type = ef;
  static constexpr size_t digits = Power;

  ef_element() = default;
  constexpr explicit ef_element(const storage_type &v) : value(v) {
    if (value & ~n)
      throw std::runtime_error("Value is not an element of the field.");
  }
  template <typename Mp>
  constexpr explicit ef_element(const ef_element<2, 1, Mp> &e)
      : value(e ? 1 : 0) {}

  static constexpr ef_element from_power(unsigned power) {
    return exp.at(power % size);
  }

  unsigned power() const { return log.at(value); }

  ef_element operator+(const ef_element &rhs) const {
    return ef_element(value ^ rhs.value);
  }

  ef_element operator-(const ef_element &rhs) const { return *this + rhs; }

  ef_element operator*(const ef_element &rhs) const {
    if (value == 0 || rhs.value == 0)
      return ef_element(0);
    return ef_element(exp.at(power() + rhs.power()));
  }

  ef_element operator/(const ef_element &rhs) const {
    if (rhs.value == 0)
      throw std::overflow_error("Divide by zero exception");

    if (value == 0)
      return ef_element(0);

    return exp.at(static_cast<size_t>(power()) -
                  static_cast<size_t>(rhs.power()) + size - 1);
  }
  ef_element inverse() const { return ef_element(1) / *this; }

  ef_element &operator*=(const ef_element &rhs) {
    *this = *this * rhs;
    return *this;
  }

  ef_element &operator+=(const ef_element &rhs) {
    *this = *this + rhs;
    return *this;
  }

  ef_element &operator++() {
    *this = *this + ef_element(1);
    /*
    if (value == 0) {
      value = 1;
    } else {
      value = exp.at(power() + 1);
    }
    */
    return *this;
  }

  ef_element operator++(int) {
    ef_element result(*this);
    ++(*this);
    return result;
  }

  bool operator<(const ef_element &rhs) const {
    /* x < 0 is false */
    if (!rhs.value)
      return false;
    /* 0 < !0 is true */
    if (!value)
      return true;
    /* sort by power */
    return power() < rhs.power();
  }
  bool operator!=(const ef_element &rhs) const { return value != rhs.value; }
  bool operator==(const ef_element &rhs) const { return value == rhs.value; }

  explicit operator bool() const { return value != 0; }
  explicit operator storage_type() const { return value; }
  explicit operator unsigned() const { return value; }
  explicit operator unsigned long() const { return value; }
  explicit operator unsigned long long() const { return value; }
  explicit operator int() const { return value; }
  explicit operator float() const { return value; }

  friend std::ostream &operator<<(std::ostream &os, const ef_element &e) {
    if (e.value == 0)
      return os << 0;
    else
      return os << "α^" << e.power();
  }
};

template <typename GF>
static constexpr std::pair<typename GF::Exp_table_type,
                           typename GF::Log_table_type>
init_tables(const long dimension, const unsigned long size,
            const uint16_t modular_polynomial) {
  using storage_type = typename GF::storage_type;
  using element_type = typename GF::element_type;

  typename GF::Exp_table_type exp;
  typename GF::Log_table_type log;

  storage_type polynomial = 1;
  for (storage_type power = 0; power < size - 1; power++) {
    log.at(polynomial) = power;
    log.at(polynomial + size) = power;

    exp.at(power) = element_type(polynomial);
    exp.at(power + size - 1) = element_type(polynomial);

    bool carry = polynomial & (1 << (dimension - 1));
    polynomial <<= 1;
    if (carry)
      polynomial ^= modular_polynomial;
  }

  log.at(0) = 0;
  log.at(size) = 0;

  exp.at(size - 1) = element_type(1);
  exp.at(2 *size - 2) = element_type(1);

  return std::make_pair(exp, log);
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
template <long Power, typename Mp>
const typename ef_element<2, Power, Mp>::Log_table_type
ef_element<2, Power, Mp>::log =
    init_tables<ef_element<2, Power, Mp> >(Power,
                                           ef_element<2, Power, Mp>::size,
                                           Mp::value).second;

template <long Power, typename Mp>
const typename ef_element<2, Power, Mp>::Exp_table_type
ef_element<2, Power, Mp>::exp =
    init_tables<ef_element<2, Power, Mp> >(Power,
                                           ef_element<2, Power, Mp>::size,
                                           Mp::value).first;
#pragma clang diagnostic pop
}

