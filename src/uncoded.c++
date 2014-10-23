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
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;

  int k;
  int dmin;
  typename std::mt19937_64::result_type seed = 0;
  // std::chrono::high_resolution_clock::now().time_since_epoch().count();

  while (1) {
    static struct option options[] = {
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
                             std::get<1>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));


  auto num_samples = [=](const double ber) {
    return std::min(1000 * 1 / ber, 10e6);
    // return std::min(base_trials * pow(10, eb_no / 2), 10e6);
  };

  auto run_algorithm = [&](const auto &algorithm) {
    const auto &name = algorithm.first;
    const auto &decoder = algorithm.second;

    std::ostringstream os;
    os << name << "_" << std::get<1>(parameters) << "_"
       << std::get<1>(parameters) << "_" << std::get<2>(parameters);
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

  const auto &uncoded = get_algorithms<0>(code).back();

  std::cout << "Running algorithm " << uncoded.first << std::endl;
  run_algorithm(uncoded);
}

