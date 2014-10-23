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
  std::cout << "--l <num>        "
            << "  "
            << "Choose code length l" << std::endl;
  std::cout << "--seed <num>     "
            << "  "
            << "Set seed of the random number generator." << std::endl;
  exit(-1);
}

int main(int argc, char *const argv[]) {
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;

  int l;
  typename std::mt19937_64::result_type seed = 0;
  // std::chrono::high_resolution_clock::now().time_since_epoch().count();

  while (1) {
    static struct option options[] = {
      { "l", required_argument, nullptr, 'l' },
      { "seed", required_argument, nullptr, 's' },
      { nullptr, 0, nullptr, 0 },
    };

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "", options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'l':
      l = strtoul(optarg, nullptr, 0);
      break;
    case 's':
      seed = strtoull(optarg, nullptr, 0);
      break;
    default:
      std::cerr << "Unkown argument: " << c << " " << std::endl;
      usage();
    }
  }

  if (!l) {
    std::cerr << "l not set." << std::endl;
    usage();
  }

  std::mt19937_64 generator;

  auto simulator =
      std::bind(evaluate, generator, std::placeholders::_1,
                std::placeholders::_2, std::placeholders::_3, l, l, 0);

  auto num_samples = [=](const double ber) {
    return std::min(10000 * 1 / ber, 10e6);
  };

  auto run_algorithm = [&](const auto &algorithm) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;

    std::ostringstream os;
    os << name << "_" << l << "_" << l << "_" << 0 << ".dat";
    std::string fname(os.str());

    if (file_exists(fname)) {
      std::cout << "File " << fname << " already exists." << std::endl;
      return;
    }

    std::ofstream file(fname.c_str(), std::ofstream::out);
    file << "eb_no " << eval_object::header() << std::endl;

    generator.seed(seed);

    double wer = 0.5;
    for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
      auto samples = num_samples(wer);
      std::cout << "Calculating for E_b/N_0 = " << eb_no << " with " << samples
                << " … ";
      std::cout.flush();
      auto start = std::chrono::high_resolution_clock::now();

      auto result = simulator(decoder, samples, eb_no);

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      std::cout << seconds << " s" << std::endl;

      wer = result.wer();

      file << std::setprecision(1) << eb_no << " ";
      file << result;
      file << std::endl;
    }
  };

  std::cout << "Running algorithm uncoded" << std::endl;
  run_algorithm(std::make_pair(std::string("uncoded"), uncoded));
}

