#pragma once

#include <algorithm>

#include "codes.h"
#include "gf/polynomial.h"
#include "math/matrix.h"

#include "hard_decision.h"
#include "soft_decision.h"

namespace cyclic {

struct coding_tag {};

struct multiplication_tag : coding_tag {};
struct division_tag : coding_tag {};
struct generator_matrix_tag : coding_tag {};

struct error_values_tag {};
struct naive_tag : error_values_tag {};
struct forney_tag : error_values_tag {};

template <typename Polynomial>
Polynomial encode(const Polynomial &g, const Polynomial &a,
                  multiplication_tag) {
  return a * g;
}

template <typename Polynomial>
Polynomial encode(const Polynomial &g, const Polynomial &a, division_tag) {
  /* a * x^k */
  auto x_k = a * Polynomial::n(g.degree());
  return x_k + (x_k % g);
}

template <typename Polynomial>
Polynomial decode(const Polynomial &g, const Polynomial &b,
                  multiplication_tag) {
  return b / g;
}

template <typename Polynomial>
Polynomial decode(const Polynomial &g, const Polynomial &b, division_tag) {
  return b / Polynomial::n(g.degree());
}

template <typename Polynomial,
          typename Element = typename Polynomial::element_type>
std::vector<Element> calculate_syndromes(const Polynomial &b,
                                         const std::vector<Element> roots) {
  std::vector<Element> syndromes;
  for (const auto &root : roots) {
    syndromes.push_back(b(root));
  }

  return syndromes;
}

/* a bit of a misnomer. should be bch_base or something. */
/* TODO: add N for shortened codes. */
template <unsigned q, typename Capability,
          typename Algorithm = peterson_gorenstein_zierler_tag,
          unsigned N = (q << 1) - 1, typename Coding = division_tag,
          typename Error = naive_tag>
class cyclic {
  static_assert(std::is_base_of<coding_tag, Coding>::value,
                "Coding must be division_tag, multiplication_tag or "
                "generator_matrix_tag");

  static_assert(
      std::is_base_of<algorithm_tag, Algorithm>::value,
      "Algorithm must be peterson_tag, berlekamp_massey_tag, or euklid_tag");

  static_assert(std::is_base_of<error_values_tag, Error>::value,
                "Error must be naive_tag or forney_tag.");

  static_assert(std::is_same<errors<Capability::value>, Capability>::value ||
                    std::is_same<dmin<Capability::value>, Capability>::value,
                "Capability has to be of type errors<> or dmin<>.");

public:
  using Galois_Field =
      gf::gf<q, typename gf::default_modular_polynomial<q>::type>;
  using Element = typename Galois_Field::element_t;
  using Polynomial = gf::polynomial<Galois_Field>;
  static constexpr unsigned n = N;
  static constexpr unsigned t = correction_capability<Capability>::value;

  static const Polynomial f;

protected:
  Polynomial g;
  Polynomial h;
  std::vector<Element> roots;

  unsigned k;
  unsigned l;
  unsigned dmin;

public:
  const double rate;

private:
  static Polynomial init_f() {
    /* x^n + 1 */
    return Polynomial::n(n) + Element(1);
  }

  /* TODO: implement & benchmark chien search */
  std::vector<Element> zeroes(const Polynomial &sigma) const {
    auto zeroes = sigma.zeroes();

    std::sort(std::begin(zeroes), std::end(zeroes));
    auto last = std::unique(std::begin(zeroes), std::end(zeroes));
    zeroes.erase(last, std::end(zeroes));

    if (zeroes.size() != sigma.degree()) {
      std::ostringstream os;
      os << "Σ(x) has to have " << sigma.degree() << " zeroes, but it has "
         << zeroes.size() << "." << std::endl;
      throw decoding_failure(os.str());
    }

    if (std::none_of(std::begin(zeroes), std::end(zeroes),
                     [](const Element &zero) { return bool(zero); }))
      throw decoding_failure("0 is zero in Σ(x)");

    return zeroes;
  }

  std::vector<unsigned>
  error_positions(const std::vector<Element> &zeroes) const {
    std::vector<unsigned> error_positions;
    for (const auto &zero : zeroes)
      error_positions.push_back(zero.power());

    return error_positions;
  }

  /* narf. don't need this for binary bch codes - we should inject a functor via
   * ctor*/
  virtual std::vector<Element> error_values(const std::vector<Element> &,
                                            const std::vector<Element> &,
                                            naive_tag) const = 0;

  virtual std::vector<Element> error_values(const std::vector<Element> &,
                                            const std::vector<Element> &,
                                            forney_tag) const = 0;

  /* TODO: select functor on signedness */
  /* is signed */
  template <typename InputSequence>
  Polynomial sequence_to_polynomial(const InputSequence &b,
                                    std::true_type) const {
    Polynomial b_;
    b_.reserve(n);
    std::transform(std::cbegin(b), std::cend(b), std::back_inserter(b_),
                   [](const auto &e) { return Element(e < 0); });
    return b_;
  }

  /* unsigned */
  template <typename InputSequence>
  Polynomial sequence_to_polynomial(const InputSequence &b,
                                    std::false_type) const {
    Polynomial b_;
    b_.reserve(n);
    std::transform(std::cbegin(b), std::cend(b), std::back_inserter(b_),
                   [](const auto &e) { return Element(e); });
    return b_;
  }

  template <typename Return_type = typename Galois_Field::storage_t,
            typename InputSequence>
  Polynomial correct_(const InputSequence &b,
                      const std::vector<unsigned> &erasures,
                      hard_decision_tag) const {
    if (b.size() != n) {
      std::ostringstream os;
      os << "Channel code word has the wrong size (" << b.size()
         << "). Expected " << n;
      throw std::runtime_error(os.str());
    }

    Polynomial b_(sequence_to_polynomial(
        b,
        typename std::is_signed<typename InputSequence::value_type>::type()));

    /* add error correction */
    auto syndromes = calculate_syndromes(b_, roots);
    bool error = std::any_of(std::begin(syndromes), std::end(syndromes),
                             [&](const auto &e) { return bool(e); });

    if (error) {
      const auto sigma_ = error_locator_polynomial<Polynomial>(
          syndromes, erasures, Algorithm());
      const auto zeroes_ = zeroes(sigma_);
      /* TODO: let this be a functor supplied by the derived class */
      const auto positions = error_positions(zeroes_);
      const auto values = error_values(syndromes, zeroes_, Error());

      auto value = std::begin(values);
      for (const auto &position : positions) {
        b_.at(position) += *value++;
      }

      syndromes = calculate_syndromes(b_, roots);

      /* declare decoding failure - Avoid decoder malfunction */
      if (std::any_of(std::begin(syndromes), std::end(syndromes),
                      [&](const auto &e) { return bool(e); }))
        throw decoding_failure("Corrected word is not a codeword");
    }

    return b_;
  }

  template <typename InputSequence, typename Tag,
            typename std::enable_if<std::is_base_of<
                soft_decision_tag, Tag>::value>::type * = nullptr>
  Polynomial correct_(const InputSequence &b,
                      const std::vector<unsigned> &erasures, Tag) const {
    using Result_type = typename Galois_Field::storage_t;
    auto copy(b);

    for (const auto &erasure : erasures)
      copy.at(erasure) = typename InputSequence::value_type(0);

    std::vector<Result_type> result =
        std::get<0>(min_sum<float, Result_type>(H_alt<int>(), copy, Tag()));
    return Polynomial(result);
  }

  static unsigned consecutive_zeroes(const Polynomial &g) {
    unsigned zeroes = 0;
    for (auto it = std::cbegin(Galois_Field()) + 1;
         it != std::cend(Galois_Field()); ++it) {
      if (g(*it))
        break;
      zeroes++;
    }
    return zeroes;
  }

public:
  cyclic(Polynomial g, std::vector<Element> roots)
      : g(g), h(f / g), roots(roots), k(g.degree()), l(n - k),
        dmin(consecutive_zeroes(g) + 1), rate((double)l / n) {
    if (dmin > n) {
      std::cout << "dmin error: " << dmin << std::endl;
      throw std::runtime_error("dmin > n");
    }

    std::cout << "f(x) = " << f << std::endl;
    std::cout << "g(x) = " << g << std::endl;
    std::cout << "h(x) = " << h << std::endl;
  }
  virtual ~cyclic() = default;
  cyclic(const cyclic &) = default;
  cyclic(cyclic &&) = default;
  cyclic &operator=(const cyclic &) = default;
  cyclic &operator=(cyclic &&) = default;

  std::string to_string() const {
    std::ostringstream os;
    os << "(" << n << ", " << l << ", " << dmin << ")-"
       << Algorithm::to_string();
    return os.str();
  }

  template <typename InputSequence, typename OutputIterator>
  void encode(const InputSequence &a, OutputIterator &&out) const {
    if (a.size() != l) {
      std::ostringstream os;
      os << "Source code word has wrong length (" << a.size() << "). Expected "
         << l;
      throw std::runtime_error(os.str());
    }

    Polynomial a_;
    a_.reserve(l);
    std::transform(std::cbegin(a), std::cend(a), std::back_inserter(a_),
                   [](const auto &e) { return Element(e); });

    auto enc = ::cyclic::encode(g, a_, Coding());

    std::transform(std::cbegin(enc), std::cend(enc), out, [](const Element &e) {
      return typename InputSequence::value_type(e);
    });
    std::fill_n(out, n - enc.size(), typename InputSequence::value_type(0));
  }

  template <typename InputSequence,
            typename Return_type = typename InputSequence::value_type>
  std::vector<Return_type> decode(const InputSequence &b,
                                  const std::vector<unsigned> &erasures =
                                      std::vector<unsigned>()) const {
    auto b_ = ::cyclic::decode(g, correct_(b, erasures, Algorithm()), Coding());

    std::vector<Return_type> r;
    r.reserve(l);
    std::transform(std::cbegin(b_), std::cend(b_), std::back_inserter(r),
                   [](const auto &e) { return Return_type(e); });
    std::fill_n(std::back_inserter(r), l - b_.size(), Return_type(0));
    return r;
  }

  /* TODO let Return_type be Galois_Field::element_t for RS codes, and anything
   * for (binary) BCH codes */
  template <typename Return_type = typename Galois_Field::storage_t,
            typename InputSequence>
  std::vector<Return_type> correct(const InputSequence &b,
                                   const std::vector<unsigned> &erasures =
                                       std::vector<unsigned>()) const {
    auto b_ = correct_(b, erasures, Algorithm());

    //std::cout << b_ << std::endl;
    std::vector<Return_type> r;
    r.reserve(n);
    std::transform(std::cbegin(b_), std::cend(b_), std::back_inserter(r),
                   [](const auto &e) { return Return_type(e); });
    std::fill_n(std::back_inserter(r), n - b_.size(), Return_type(0));
    return r;
  }

  template <typename T> matrix<T> H() const {
    std::vector<T> row(n, T(0));
    std::transform(std::crbegin(h), std::crend(h), std::begin(row),
                   [](const Element &e) { return T(e); });

    matrix<T> control(row);

    std::generate_n(std::back_inserter(control), k - 1, [&]() {
      std::rotate(std::rbegin(row), std::rbegin(row) + 1, std::rend(row));
      return row;
    });

    return control;
  }

  template <typename T> matrix<T> H_alt() const {
    std::vector<Element> row;
    for (unsigned power = 0; power < n; power++)
      row.push_back(Element::from_power(0));
    matrix<Element> tmp(t, n);

    for (unsigned row = 0; row < t; row++) {
      for (unsigned col = 0; col < n; col++) {
        const unsigned row_power = 2 * row + 1;
        tmp.at(row).at(col) = Element::from_power(col * row_power);
      }
    }

    matrix<T> control(0, n);

    for (const auto &row : tmp) {
      for (size_t mask = 1, digits = Element::digits; digits;
           mask <<= 1, digits--) {
        std::vector<T> r;
        r.reserve(n);
        for (const auto &element : row)
          r.push_back(bool(static_cast<size_t>(element) & mask));
        control.push_back(r);
      }
    }

    return control;
  }
};

template <unsigned q, typename Capability, typename Algorithm, unsigned N,
          typename Coding, typename Error>
const typename cyclic<q, Capability, Algorithm, N, Coding, Error>::Polynomial
cyclic<q, Capability, Algorithm, N, Coding, Error>::f =
    cyclic<q, Capability, Algorithm, N, Coding, Error>::init_f();
}
