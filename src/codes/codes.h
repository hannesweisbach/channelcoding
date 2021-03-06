#pragma once

#include <stdexcept>
#include <numeric>
#include <algorithm>

template <unsigned e> struct errors {
  static constexpr unsigned value = e;
};

template <unsigned d> struct dmin {
  static constexpr unsigned value = d;
};

template <typename T> struct correction_capability {
  /* Be always false, but 'false' will not compile. */
  static_assert(sizeof(T) == 0, "Unsupported type to specify error correction "
                                "capabilities. Use errors<> or dmin<>.");
};
template <unsigned v> struct correction_capability<dmin<v> > {
  static constexpr unsigned value = (v - 1) / 2;
};

template <unsigned v> struct correction_capability<errors<v> > {
  static constexpr unsigned value = v;
};

class decoding_failure : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
  virtual ~decoding_failure();
  decoding_failure(const decoding_failure &);
  decoding_failure(decoding_failure &&);
  decoding_failure &operator=(const decoding_failure &);
  decoding_failure &operator=(decoding_failure &&);
};

struct algorithm_tag {};

struct hard_decision_tag : algorithm_tag {};
struct soft_decision_tag : algorithm_tag {};

template <typename T, typename InputIterator, typename OutputIterator>
void hard_decision(InputIterator first, InputIterator last,
                   OutputIterator out) {
#if 0
  // TODO(hannes): spezialize std::is_signed<> for ef_element.
  static_assert(!std::is_signed<InputIterator::value_type>::value,
                "U has to be unsigned");
#endif
  std::transform(first, last, out, [](const auto &e) { return T(e < 0); });
}
