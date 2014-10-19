#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <cstdlib>

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
      bit_errors += length;
    }
  }
  /* word error rate */
  std::ostringstream os;
  os << std::scientific;
  // os << std::setprecision(4) << eb_n0 << " ";
  // os << std::setprecision(9) << fk_corr / (float)samples << " ";
  os << std::setprecision(12) << bit_errors / (double)(samples * length) << " ";
  // os << std::setprecision(12) << reconstruction_failures / (double)samples;

  return os.str();
}

int main(int argc, const char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr unsigned base_trials = 10000;

  auto num_samples = [=](const double eb_no) {
    return std::min(base_trials * pow(10, eb_no / 2), 10e6);
  };

  const float alpha_start = (argc > 1) ? strtof(argv[1], nullptr) : 0.1f;
  const float alpha_max = (argc > 2) ? strtof(argv[2], nullptr) : 1.0f;
  const float alpha_step = (argc > 3) ? strtof(argv[3], nullptr) : 0.01f;

  const float beta_start = (argc > 4) ? strtof(argv[4], nullptr) : 0.1f;
  const float beta_max = (argc > 5) ? strtof(argv[5], nullptr) : 1.0f;
  const float beta_step = (argc > 6) ? strtof(argv[6], nullptr) : 0.01f;

  const float eb_no = (argc > 7) ? strtof(argv[7], nullptr) : 0.0f;

  std::cout << "Calculating for eb_no: " << eb_no << std::endl;
  std::cout << "alpha: " << alpha_start << " to " << alpha_max << " in "
            << alpha_step << " steps" << std::endl;
  std::cout << "beta:  " << beta_start << " to " << beta_max << " in "
            << beta_step << " steps" << std::endl;

  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  bch code(5, 0x25, 7);
  // bch code(6, 0x45, 7);

  const char *fname;
#if 0
  char fbuf[] = "logs/optimize.XXXXXX";
  mktemp(fbuf);
  fname = fbuf;
#else
  std::ostringstream os;
  os << "logs/ebno_" << eb_no;
  fname = os.str().c_str();
#endif

  std::ofstream file(fname, std::ofstream::out);
  std::cout << "Writing to file " << fname << std::endl;
  file << "#eb_no " << eb_no << std::endl;
  file << "#alpha: " << alpha_start << " " << alpha_max << " " << alpha_step
       << std::endl;
  file << "#beta : " << beta_start << " " << beta_max << " " << beta_step
       << std::endl;

  file << std::scientific;

  for (float alpha = alpha_start; alpha < alpha_max; alpha += alpha_step) {
    for (float beta = beta_start; beta < beta_max; beta += beta_step) {
      const size_t N = num_samples(eb_no);

      std::cout << "alpha = " << alpha << " beta = " << beta << " N = " << N;
      std::cout.flush();

      auto start = std::chrono::high_resolution_clock::now();
      auto f =
          std::bind(nms_2d<max_iterations, float, float>, std::ref(code.H()),
                    std::placeholders::_1, alpha, beta);
      file << generate_overview(generator, f, N, eb_no) << " ";

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

      std::cout << seconds << " s" << std::endl;
    }

    file << std::endl;
  }
}
