#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>

#include "eval.h"
#include "util.h"

std::string generate_overview(std::mt19937 &generator, const decoder_t &decoder,
                              const size_t samples, const float eb_n0,
                              const size_t n, const size_t l, const size_t fk) {
  const float R = (float)l / n;

  std::vector<float> b(l);

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
      bit_errors += l;
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
  os << std::setprecision(12) << bit_errors / (samples_ * l);

  return os.str();
}

int main(int argc, const char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;
  constexpr unsigned base_trials = 10000;

  bch code(5, 0x25, 7);
  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  auto parameters = code.parameters();
  auto simulator = std::bind(generate_overview, generator,
                             std::placeholders::_1, std::placeholders::_2,
                             std::placeholders::_3, std::get<0>(parameters),
                             std::get<1>(parameters), std::get<2>(parameters));

  auto algorithms = get_algorithms<max_iterations>(code);

  auto num_samples = [=](const double eb_no) {
    return std::min(base_trials * pow(10, eb_no / 2), 10e6);
  };

  auto run_algorithm = [&](const auto &algorithm) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;

    if (file_exists(name)) {
      std::cout << "File " << name << " already exists." << std::endl;
      return;
    }

    std::ofstream file(name.c_str(), std::ofstream::out);
    file << "eb_no "
         << "reconstruction_failures "
         << "reconstruction_errors "
         << "wer "
         << "fk_rate "
         << "ber" << std::endl;
    for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
      auto samples = num_samples(eb_no);
      std::cout << "Calculating for E_b/N_0 = " << eb_no << " with " << samples
                << " … ";
      std::cout.flush();
      auto start = std::chrono::high_resolution_clock::now();

      file << eb_no << " ";
      file << simulator(decoder, samples, eb_no);
      file << std::endl;

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

      std::cout << seconds << " s" << std::endl;
    }
  };

  if (argc < 2) {
    std::cout << "Running all algorithms" << std::endl;
    for (const auto &algorithm : algorithms) {
      run_algorithm(algorithm);
    }
  } else {
    for (int arg = 1; arg < argc; arg++) {
      unsigned long index = strtoul(argv[arg], nullptr, 0);
      if (index < algorithms.size()) {
        std::cout << "Running algorithm " << algorithms.at(index).first
                  << std::endl;
        run_algorithm(algorithms.at(index));
      } else {
        unsigned index = 0;
        std::cout << "Index " << index
                  << " is out of range. Choose from:" << std::endl;
        for (const auto &algorithm : algorithms)
          std::cout << "[" << index++ << "] " << algorithm.first << std::endl;
      }
    }
  }
}
