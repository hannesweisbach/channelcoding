#include "eval.h"

#include <cmath>
#include <vector>

uint64_t primitive_polynomial(unsigned degree) {
  static const std::vector<uint64_t> primitive_polynomials(
      { 0, 0x3, 0x7, 0xb, 0x13, 0x25, 0x43 });
  return primitive_polynomials.at(degree);
}

static double sigma(const double eb_n0, const double R) {
  return 1.0f / sqrt((2 * R * pow(10, eb_n0 / 10.0)));
}

std::tuple<std::vector<int>, std::vector<float>, unsigned>
pzg_wrapper(const bch &code, const std::vector<float> &b) {
  std::vector<int> hard;
  hard.reserve(b.size());
  std::transform(std::cbegin(b), std::cend(b), std::back_inserter(hard),
                 [](const auto &bit) { return bit < 0; });
  auto b_corr = code.correct_peterson(hard);
  return std::make_tuple(b_corr, std::vector<float>(), 0);
}

std::tuple<std::vector<int>, std::vector<float>, unsigned>
bm_wrapper(const bch &code, const std::vector<float> &b) {
  std::vector<int> hard;
  hard.reserve(b.size());
  std::transform(std::cbegin(b), std::cend(b), std::back_inserter(hard),
                 [](const auto &bit) { return bit < 0; });
  auto b_corr = code.correct_bm(hard);
  return std::make_tuple(b_corr, std::vector<float>(), 0);
}

std::ostream &operator<<(std::ostream &os, const eval_object &e) {
  os << std::scientific;
  const double samples_ = e.samples;
  os << std::setprecision(12) << e.reconstruction_failures / samples_ << " ";
  os << std::setprecision(12) << e.reconstruction_errors / samples_ << " ";
  const auto wrong_words = e.reconstruction_failures + e.reconstruction_errors;
  os << std::setprecision(12) << wrong_words / samples_ << " ";
  os << std::setprecision(12) << e.fk_corr / samples_ << " ";
  os << std::setprecision(12) << e.bit_errors / (samples_ * e.n);
  return os;
}

double eval_object::ber() const { return bit_errors / (double)(samples * n); }

eval_object evaluate(std::mt19937_64 &generator, const decoder_t &decoder,
                     const size_t samples, const float eb_n0, const unsigned n,
                     const unsigned l, const unsigned fk) {
  const float R = (float)l / n;

  unsigned reconstruction_failures = 0;
  unsigned reconstruction_errors = 0;
  unsigned fk_corr = 0;
  unsigned bit_errors = 0;

  std::vector<float> b(n);
  std::normal_distribution<float> random(1.0, sigma(eb_n0, R));
  auto noise_gen = std::bind(std::ref(random), std::ref(generator));

  for (size_t sample = 0; sample < samples; sample++) {
    std::generate(std::begin(b), std::end(b), noise_gen);
    try {
      auto result = decoder(b);
      const auto &b_corr = std::get<0>(result);
      auto wrong_bits = std::count_if(std::cbegin(b_corr), std::cend(b_corr),
                                      [](const auto &bit) { return bit != 0; });
      if (wrong_bits) {
        reconstruction_errors++;
        bit_errors += wrong_bits;
      } else {
        auto d = std::count_if(std::cbegin(b), std::cend(b),
                               [](const auto &bit) { return bit < 0; });
        if (d > fk)
          fk_corr++;
      }
    }
    catch (const decoding_failure &) {
      reconstruction_failures++;
      bit_errors += n;
    }
  }

  return { eb_n0,      reconstruction_failures, reconstruction_errors, fk_corr,
           bit_errors, samples,                 n };
}

