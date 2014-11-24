#pragma once

#include <string>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <sstream>

class uncoded {
  template <typename Out, typename In,
            typename std::enable_if<std::is_signed<In>::value>::type * =
                nullptr>
  static Out binary_output(const In &e) {
    return Out(e < 0);
  }

  template <typename Out, typename In,
            typename std::enable_if<!std::is_signed<In>::value>::type * =
                nullptr>
  static Out binary_output(const In &e) {
    return Out(e);
  }

public:
  static constexpr double rate = 0.5;
  const unsigned n;

  explicit uncoded(const unsigned l) : n(l) {}

  std::string to_string() const {
    std::ostringstream os;
    os << n << "-uncoded";
    return os.str();
  }

  template <typename Return_type = uint8_t, typename InputSequence>
  std::vector<Return_type> correct(const InputSequence &b) const {
    std::vector<Return_type> r;
    r.reserve(n);
    std::transform(
        std::cbegin(b), std::cend(b), std::back_inserter(r),
        &binary_output<Return_type, typename InputSequence::value_type>);
    std::fill_n(std::back_inserter(r), n - b.size(), Return_type(0));
    return r;
  }
};

constexpr double uncoded::rate;
