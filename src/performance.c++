#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>

#include "bch.h"
#include "iterative.h"
#include "util.h"

static double sigma(const double eb_n0, const double R) {
  return sqrt(1.0f / (2 * R * pow(10, eb_n0 / 10.0)));
}

using decoder_t =
    std::function<std::tuple<std::vector<int>, std::vector<float>, unsigned>(
        std::vector<float>)>;

std::string generate_overview(std::mt19937 &generator, const decoder_t &decoder,
                              const size_t samples, const float eb_n0) {
  constexpr int fk = 3;
  constexpr size_t length = 31;
  constexpr float R = 16.0 / length;

  std::vector<float> b(length);

  std::normal_distribution<float> random(1.0, sigma(eb_n0, R));
  auto noise_gen = std::bind(std::ref(random), std::ref(generator));

  unsigned reconstruction_failures = 0;
  unsigned reconstruction_errors = 0;
  unsigned fk_corr = 0;
  unsigned bit_errors = 0;

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
    catch (...) {
      reconstruction_failures++;
    }
  }
  /* word error rate */
  std::ostringstream os;
  os << std::scientific;
  const double samples_ = samples;
  os << std::setprecision(12) << reconstruction_failures / samples_ << " ";
  os << std::setprecision(12) << reconstruction_errors / samples_ << " ";
  const auto wrong_words = reconstruction_failures + reconstruction_errors;
  os << std::setprecision(12) << wrong_words / samples_ << " ";
  os << std::setprecision(12) << fk_corr / samples_ << " ";
  os << std::setprecision(12) << bit_errors / (samples_ * length);

  return os.str();
}

std::tuple<std::vector<int>, std::vector<float>, unsigned>
    pzg_wrapper(const bch &code, const std::vector<float> &b) {
  std::vector<int> hard;
  hard.reserve(b.size());
  std::transform(std::cbegin(b), std::cend(b), std::back_inserter(hard),
                 [](const auto &bit) { return -2 * bit + 1; });
  auto b_corr = code.correct_peterson(hard);
  return std::make_tuple(b_corr, std::vector<float>(), 0);
}

int main() {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;
  constexpr unsigned base_trials = 10000;

  bch code(5, 0x25, 7);
  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  auto num_samples = [=](const double eb_no) {
    return base_trials * pow(10, eb_no / 2);
  };

  std::vector<std::pair<std::string, decoder_t> > algorithms;
  algorithms.emplace_back("ms", std::bind(min_sum<max_iterations, float, float>,
                                          std::cref(code.H()),
                                          std::placeholders::_1));
  algorithms.emplace_back(
      "scms1", std::bind(scms1<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));
  algorithms.emplace_back(
      "scms2", std::bind(scms2<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));
  algorithms.emplace_back(
      "pzg", std::bind(pzg_wrapper, std::cref(code), std::placeholders::_1));

  for (const auto &algorithm : algorithms) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;
    std::ofstream file(name.c_str(), std::ofstream::out);
    file << "eb_no "
         << "reconstruction_failures "
         << "reconstruction_errors "
         << "wer "
         << "fk_rate "
         << "ber" << std::endl;
    for (double eb_no = 0; eb_no < 10; eb_no += eb_no_step) {
      file << eb_no << " ";
      file << generate_overview(generator, decoder, num_samples(eb_no), eb_no)
           << " ";
      file << std::endl;
    }
  }
}
