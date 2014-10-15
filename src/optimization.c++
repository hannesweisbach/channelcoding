#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>

#include "bch.h"
#include "iterative.h"
#include "util.h"

static double sigma(const double eb_n0, const double R) {
  return sqrt(1.0f / (2 * R * pow(10, eb_n0 / 10.0)));
}

using decoder_t =
    std::function<std::tuple<std::vector<int>, std::vector<float>, unsigned>(
        std::vector<float>)>;

double generate_overview(std::mt19937 &generator, decoder_t decoder,
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
  // os << std::setprecision(4) << eb_n0 << " ";
  // os << std::setprecision(9) << failures / (float)samples << " ";
  // os << std::setprecision(9) << fk_corr / (float)samples << " ";
  // os << std::setprecision(9) << bit_errors / (float)(samples * length);

  return bit_errors / (double)(samples * length);
}

int main() {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_step = 0.1f;
  constexpr float alpha_step = 0.01f;
  constexpr unsigned N = 10000;

  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  bch code(5, 0x25, 7);

  std::ofstream file("test.output", std::ofstream::out);

  /* header */
  file << "eb_no ";
  for (float alpha = 0; alpha < 1; alpha += alpha_step) {
    std::ostringstream col_name;
    col_name << "alpha_" << alpha;
    file << col_name.str() << " ";
  }
  file << std::endl;

  file << std::scientific;

  for (double eb_no = 0; eb_no < 10; eb_no += eb_no_step) {
    file << eb_no << " ";
    auto f = std::bind(min_sum<max_iterations, float, float>, std::cref(code.H()),
                       std::placeholders::_1);
    file << generate_overview(generator, f, N, eb_no) << " ";
    file << std::endl;
  }

  for (double eb_no = 0; eb_no < 10; eb_no += eb_no_step) {
    file << eb_no << " ";
    for (float alpha = 0; alpha < 1; alpha += alpha_step) {
      auto f = std::bind(nms<max_iterations, float, float>, std::ref(code.H()),
                         std::placeholders::_1, alpha);
      file << generate_overview(generator, f, N, eb_no) << " ";
      // std::cout << eb_no << " " << alpha << std::endl;
    }
    file << std::endl;
  }
}
