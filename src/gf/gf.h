#pragma once

#include <limits>
#include <cstdint>
#include <array>
#include <iterator>

#include <iostream>
#include <iomanip>

namespace gf {

/* uint16_t allows galois fields up to 2^15 to be used. */
template <uint16_t poly> struct modular_polynomial {
  static constexpr uint16_t value = poly;
};

namespace detail {
static constexpr std::array<unsigned, 9> modular_polynomials = {
  { 0, 0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d }
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

template <unsigned q> struct default_modular_polynomial {
  static_assert(q > 0, "GF(2^0) does not make sense. Choose q > 0.");
  static_assert(q < 9, "modular polynomial for GF(q), q > 2^8 have to "
                       "specified manually (or the list of default modular "
                       "polynomials has to be extended.");
  using type = modular_polynomial<std::get<q>(detail::modular_polynomials)>;
};



/* Galois field over 2^q.
 * TODO implement Galois field over prime^q.
 */
/* Think about letting the field be an inner class of gf_element.
 * This way, using gf_element as template parameter would be cleaner.
 * Iteration over all field elements could be done by:
 * a) letting gf_element::operator++ increment the power of alpha
 * b) iterating over all powers and using gf_element::from_power
 * c) having an inner class field in gf_element with cbegin/cend
 * a) isn't nice, because the field element 0 is not a power of alpha.
 * b) isn't nice, because operator++ has weird semantics, but it would sort-of
 *    do what one would expect
 * c) seems the best solution ...
 */

template <unsigned q, typename Modular_Polynomial =
                          typename default_modular_polynomial<q>::type>
class gf {
  static_assert(
      std::is_same<Modular_Polynomial,
                   modular_polynomial<Modular_Polynomial::value> >::value,
      "Modular_Polynomial has to be of type modular_polynomial<>.");
  class gf_element;

public:
  using element_t = gf_element;
  using element_type = gf_element;
  /* TODO: think about using std::bitset. */
  using storage_t = typename detail::uint_t<q>::type;

  /* size of the field */
  static constexpr size_t size = (1 << q);
  static constexpr unsigned modular_polynomial = Modular_Polynomial::value;

private:
  static_assert(q, "q must be non-zero.");
  static_assert(modular_polynomial & (1 << q),
                "q-th bit of modular polynomial has to be set.");
  static_assert(
      (modular_polynomial & ~(size - 1)) <= size,
      "q-th bit has to be the highest set bit in the modular polynomial.");

  /* number of non-zero elements in the field */
  static constexpr size_t n = size - 1;

  using Log_table_type = std::array<storage_t, 2 * size>;
  using Exp_table_type = std::array<gf_element, 2 * size>;

  template <typename GF, unsigned>
  friend constexpr std::pair<typename GF::Exp_table_type,
                             typename GF::Log_table_type>
  init_tables();

  static Log_table_type log;
  static Exp_table_type exp;

  class gf_element {
    storage_t value;

  public:
    using storage_t = storage_t;
    using storage_type = storage_t;
    using field_type = gf;

    static constexpr size_t digits = q;
    constexpr explicit gf_element() : value(0) {}
    constexpr explicit gf_element(storage_t v) : value(v) {
      if (value & ~n)
        throw std::runtime_error("Value is not an element of the field.");
    }
    static constexpr gf_element from_power(unsigned power) {
      return exp.at(power % size);
    }

    unsigned power() const { return log.at(value); }

    gf_element operator+(const gf_element &rhs) const {
      return gf_element(value ^ rhs.value);
    }

    gf_element operator-(const gf_element &rhs) const { return *this + rhs; }

    gf_element operator*(const gf_element &rhs) const {
      if (value == 0 || rhs.value == 0)
        return gf_element(0);
      return gf_element(exp.at(power() + rhs.power()));
    }

    gf_element operator/(const gf_element &rhs) const {
      if (rhs.value == 0)
        throw std::overflow_error("Divide by zero exception");

      if (value == 0)
        return gf_element(0);

      return exp.at(static_cast<size_t>(power()) -
                    static_cast<size_t>(rhs.power()) + size - 1);
    }

    gf_element &operator*=(const gf_element &rhs) {
      *this = *this * rhs;
      return *this;
    }

    gf_element &operator+=(const gf_element &rhs) {
      *this = *this + rhs;
      return *this;
    }

    gf_element &operator++() {
      *this = *this + one;
      /*
      if (value == 0) {
        value = 1;
      } else {
        value = exp.at(power() + 1);
      }
      */
      return *this;
    }

    gf_element operator++(int) {
      gf_element result(*this);
      ++(*this);
      return result;
    }

    gf_element inverse() const { return one / *this; }
    gf_element zero() const;

    bool operator<(const gf_element &rhs) const {
      /* x < 0 is false */
      if (!rhs.value)
        return false;
      /* 0 < !0 is true */
      if (!value)
        return true;
      /* sort by power */
      return power() < rhs.power();
    }
    bool operator!=(const gf_element &rhs) const { return value != rhs.value; }
    bool operator==(const gf_element &rhs) const { return value == rhs.value; }

    explicit operator bool() const { return value != 0; }
    explicit operator storage_t() const { return value; }
    explicit operator unsigned() const { return value; }
    explicit operator unsigned long() const { return value; }
    explicit operator unsigned long long() const { return value; }
    explicit operator int() const { return value; }
    explicit operator float() const { return value; }

    friend std::ostream &operator<<(std::ostream &os, const gf_element &e) {
      if (e.value == 0)
        return os << 0;
      else
        return os << "α^" << e.power();
    }
  };

public:
  gf() = default;
  using const_iterator = typename Exp_table_type::const_iterator;
  using const_reverse_iterator =
      typename Exp_table_type::const_reverse_iterator;

  static const element_t zero;
  static const element_t one;

  static constexpr element_t element(unsigned power) { return exp.at(power); }

  const_iterator begin() const { return std::begin(exp); }
  const_iterator end() const { return std::begin(exp) + size; }
  const_iterator cbegin() const noexcept { return std::cbegin(exp); }
  const_iterator cend() const noexcept { return std::cbegin(exp) + size; }
  const_reverse_iterator rbegin() const { return std::rbegin(exp) + size; }
  const_reverse_iterator rend() const { return std::rend(exp); }
  const_reverse_iterator crbegin() const noexcept {
    return std::crbegin(exp) + size;
  }
  const_reverse_iterator crend() const noexcept { return std::crend(exp); }

  friend std::ostream &operator<<(std::ostream &os, const gf &) {
    for (int i = 0; i < log.size(); i++)
      std::cout << i << " " << static_cast<int>(log.at(i)) << std::endl;

    for (int i = 0; i < exp.size(); i++)
      std::cout << i << " "
                << static_cast<int>(static_cast<storage_t>(exp.at(i)))
                << std::endl;

    os << "  α    poly   α    poly" << std::endl;
    for (unsigned power = 0; power < size; power++) {
      auto element = exp.at(power);
      storage_t bits = static_cast<storage_t>(element);
      if (bits && (power != size - 1) && (power != log.at(bits))) {
        std::cout << "  Expected power " << power
                  << ", but got power: " << static_cast<unsigned>(log.at(bits))
                  << std::endl;
        std::cout << "  Expected element " << element << ", but got "
                  << exp.at(power) << std::endl;
      }
      os << "α^" << std::setw(3) << std::left << power << ": " << std::hex
         << static_cast<unsigned>(bits) << std::dec << "    ";
      os << "α^" << std::setw(3) << std::left << power + size << ": "
         << std::hex
         << static_cast<unsigned>(static_cast<storage_t>(exp.at(power + size)))
         << std::dec;
      if (power < size - 1)
        os << std::endl;
    }
    return os;
  }
};

template <typename GF, unsigned q>
static constexpr std::pair<typename GF::Exp_table_type,
                           typename GF::Log_table_type>
init_tables() {
  using storage_t = typename GF::storage_t;
  using element_t = typename GF::element_t;

  constexpr size_t size = GF::size;

  typename GF::Exp_table_type exp;
  typename GF::Log_table_type log;

  storage_t polynomial = 1;
  for (storage_t power = 0; power < size - 1; power++) {
    log.at(polynomial) = power;
    log.at(polynomial + size) = power;

    exp.at(power) = element_t(polynomial);
    exp.at(power + size - 1) = element_t(polynomial);

    bool carry = polynomial & (1 << (q - 1));
    polynomial <<= 1;
    if (carry)
      polynomial ^= GF::modular_polynomial;
  }

  log.at(0) = 0;
  log.at(size) = 0;

  exp.at(size - 1) = element_t(1);
  exp.at(2 * size - 2) = element_t(1);

  return std::make_tuple(exp, log);
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"

template <unsigned q, typename Modular_Polynomial>
typename gf<q, Modular_Polynomial>::Log_table_type
gf<q, Modular_Polynomial>::log =
    init_tables<gf<q, Modular_Polynomial>, q>().second;

template <unsigned q, typename Modular_Polynomial>
typename gf<q, Modular_Polynomial>::Exp_table_type
gf<q, Modular_Polynomial>::exp =
    init_tables<gf<q, Modular_Polynomial>, q>().first;

#pragma clang diagnostic pop

template <unsigned q, typename Modular_Polynomial>
const typename gf<q, Modular_Polynomial>::element_t
gf<q, Modular_Polynomial>::one =
    typename gf<q, Modular_Polynomial>::element_t(1);

template <unsigned q, typename Modular_Polynomial>
const typename gf<q, Modular_Polynomial>::element_t
gf<q, Modular_Polynomial>::zero =
    typename gf<q, Modular_Polynomial>::element_t(0);
}
