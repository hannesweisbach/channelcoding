#include "eval.h"

#include <cmath>
#include <vector>
#include <chrono>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>

uint64_t primitive_polynomial(unsigned degree) {
  static const std::vector<uint64_t> primitive_polynomials(
      { 0, 0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83 });
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

std::tuple<std::vector<int>, std::vector<float>, unsigned>
uncoded(const std::vector<float> &b) {
  std::vector<int> hard;
  hard.reserve(b.size());
  std::transform(std::cbegin(b), std::cend(b), std::back_inserter(hard),
                 [](const auto &bit) { return bit < 0; });
  return std::make_tuple(hard, std::vector<float>(), 0);
}

std::string eval_object::header() {
  return "reconstruction_failures "
         "reconstruction_errors "
         "wer "
         "fk_rate "
         "ber ";
}

std::ostream &operator<<(std::ostream &os, const eval_object &e) {
  os << std::scientific;
  const double samples_ = e.samples;
  os << std::setprecision(12) << e.reconstruction_failures / samples_ << " ";
  os << std::setprecision(12) << e.reconstruction_errors / samples_ << " ";
  os << std::setprecision(12) << e.wer() << " ";
  os << std::setprecision(12) << e.fk_corr / samples_ << " ";
  os << std::setprecision(12) << e.ber() << " ";
  return os;
}

double eval_object::ber() const { return bit_errors / (double)(samples * n); }
double eval_object::wer() const {
  const auto wrong_words = reconstruction_failures + reconstruction_errors;
  return wrong_words / (double)(samples);
}

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

static void usage() {
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
  std::cout << "--seed-time      "
            << "  "
            << "Use the current time as seed for the random number generator."
            << std::endl;

  exit(-1);
}

std::tuple<std::vector<int>, int, int, uint64_t>
parse_options(const int argc, char *const argv[]) {

  std::vector<int> indices;
  int k;
  int dmin;
  uint64_t seed = 0;

  while (1) {
    static struct option options[] = {
      { "algorithm", required_argument, nullptr, 'a' },
      { "k", required_argument, nullptr, 'k' },
      { "dmin", required_argument, nullptr, 'd' },
      { "seed", required_argument, nullptr, 's' },
      { "seed-time", no_argument, nullptr, 't' },
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
    case 't':
      seed =
          std::chrono::high_resolution_clock::now().time_since_epoch().count();
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

  return std::make_tuple(indices, k, dmin, seed);
}
