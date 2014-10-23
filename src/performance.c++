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

void usage() {
  int index = 0;
  bch code(5, 0x25, 7);
  std::cout << "--algorithm <num>"
            << "  "
            << "Choose algorithm:" << std::endl;
  for (const auto &algorithm : get_algorithms<0>(code))
    std::cout << "  [" << index++ << "] " << algorithm.first << std::endl;
  std::cout << "--k <num>        "
            << "  "
            << "Choose code length n = 2^k - 1." << std::endl;
  std::cout << "--dmin <num      "
            << "  "
            << "Choose dmin of the code." << std::endl;
  std::cout << "--seed <num>     "
            << "  "
            << "Set seed of the random number generator." << std::endl;

  exit(-1);
}

int main(int argc, char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;
  constexpr unsigned base_trials = 10000;

  std::vector<int> indices;
  int k;
  int dmin;
  typename std::mt19937_64::result_type seed = 0;
  // std::chrono::high_resolution_clock::now().time_since_epoch().count();

  while (1) {
    static struct option options[] = {
      { "algorithm", required_argument, nullptr, 'a' },
      { "k", required_argument, nullptr, 'k' },
      { "dmin", required_argument, nullptr, 'd' },
      { "seed", required_argument, nullptr, 's' },
      { nullptr, 0, nullptr, 0 },
    };

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "", options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'a':
      indices.push_back(strtoul(optarg, nullptr, 0));
      break;
    case 'k':
      k = strtoul(optarg, nullptr, 0);
      break;
    case 'd':
      dmin = strtoul(optarg, nullptr, 0);
      break;
    case 's':
      seed = strtoull(optarg, nullptr, 0);
      break;
    default:
      std::cerr << "Unkown argument: " << c << " " << std::endl;
      usage();
    }
  }

  bool fail = false;

  if (indices.size() == 0) {
    std::cerr << "Algorithm not set." << std::endl;
    fail = true;
  }

  if (!k) {
    std::cerr << "k not set." << std::endl;
    fail = true;
  }

  if (!dmin) {
    std::cerr << "dmin not set." << std::endl;
    fail = true;
  }

  if (fail)
    usage();

  bch code(k, primitive_polynomial(k), dmin);
  std::mt19937_64 generator;

  auto parameters = code.parameters();
  auto simulator = std::bind(evaluate, generator, std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             std::get<0>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));

  auto algorithms = get_algorithms<max_iterations>(code);

  auto num_samples = [=](const double ber) {
    return std::min(1000 * 1 / ber, 10e6);
    // return std::min(base_trials * pow(10, eb_no / 2), 10e6);
  };

  auto run_algorithm = [&](const auto &algorithm) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;

    std::ostringstream os;
    os << name << "_" << std::get<0>(parameters) << "_"
       << std::get<1>(parameters) << "_" << std::get<2>(parameters);
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
