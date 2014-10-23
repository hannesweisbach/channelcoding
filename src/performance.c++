#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>

#include "eval.h"
#include "util.h"

int main(int argc, char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;
  constexpr unsigned base_trials = 10000;

  int k;
  int dmin;
  uint64_t seed;
  std::vector<int> indices;

  std::tie(indices, k, dmin, seed) = parse_options(argc, argv);

  bch code(k, primitive_polynomial(k), dmin);
  std::mt19937_64 generator;

  auto parameters = code.parameters();
  auto simulator = std::bind(evaluate, generator, std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             std::get<0>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));

  auto algorithms = get_algorithms<max_iterations>(code);

  auto run_algorithm = [&](const auto &algorithm) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;

    std::ostringstream os;
    os << name << "_" << std::get<0>(parameters) << "_"
       << std::get<1>(parameters) << "_" << std::get<2>(parameters) << ".dat";
    std::string fname(os.str());

    if (file_exists(fname)) {
      std::cout << "File " << fname << " already exists." << std::endl;
      return;
    }

    std::ofstream file(fname.c_str(), std::ofstream::out);
    file << "eb_no " << eval_object::header() << std::endl;

    generator.seed(seed);

    double ber = 0.5;
    for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
      auto samples = iterations(ber);
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

  for (const auto &index : indices) {
    if (index < algorithms.size()) {
      std::cout << "Running algorithm " << algorithms.at(index).first
                << std::endl;
      run_algorithm(algorithms.at(index));
    } else {
      std::cout << "Index " << index
                << " is out of range. Choose from:" << std::endl;
      unsigned i = 0;
      for (const auto &algorithm : algorithms)
        std::cout << "[" << i++ << "] " << algorithm.first << std::endl;
    }
  }
}
