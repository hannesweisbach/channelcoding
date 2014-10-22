#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>

#include "eval.h"
#include "util.h"

int main(int argc, const char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;
  constexpr unsigned base_trials = 10000;

  auto seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

  bch code(5, 0x25, 7);
  std::mt19937_64 generator;

  auto parameters = code.parameters();
  auto simulator = std::bind(evaluate, generator, std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             std::get<0>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));

  auto algorithms = get_algorithms<max_iterations>(code);

  auto num_samples = [=](const double ber) {
    return std::min(1000 * 1/ber, 10e6);
    //return std::min(base_trials * pow(10, eb_no / 2), 10e6);
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

    generator.seed(seed);

    double ber = 1;
    for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
      auto samples = num_samples(ber);
      std::cout << "Calculating for E_b/N_0 = " << eb_no << " with " << samples
                << " … ";
      std::cout.flush();
      auto start = std::chrono::high_resolution_clock::now();

      auto result = simulator(decoder, samples, eb_no);

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      std::cout << seconds << " s" << std::endl;

      ber = result.ber();

      file << std::setprecision(1) << eb_no << " ";
      file << result;
      file << std::endl;
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
