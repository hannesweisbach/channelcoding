#include <vector>
#include <tuple>
#include <iostream>
#include <utility>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>

#include "simulation.h"

#include "codes/bch.h"
#include "codes/rs.h"

std::vector<decoder> decoders {
#if 0
  cyclic::primitive_bch<5, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::berlekamp_massey_tag>(),
  
  cyclic::primitive_bch<5, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),

  cyclic::primitive_bch<5, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::euklid_tag>(),
#endif
  cyclic::primitive_bch<5, dmin<3>, min_sum_tag<50> >(),
      cyclic::primitive_bch<5, dmin<5>, min_sum_tag<50> >(),
      cyclic::primitive_bch<5, dmin<7>, min_sum_tag<50> >(),
      cyclic::primitive_bch<5, dmin<9>, min_sum_tag<50> >(),
      cyclic::primitive_bch<6, dmin<3>, min_sum_tag<50> >(),
      cyclic::primitive_bch<6, dmin<5>, min_sum_tag<50> >(),
      cyclic::primitive_bch<6, dmin<7>, min_sum_tag<50> >(),
#if 0
  cyclic::primitive_bch<6, dmin<9>, min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<3>, min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<5>, min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<7>, min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<9>, min_sum_tag<50>>(),

  cyclic::primitive_bch<5, dmin<3>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<5, dmin<5>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<5, dmin<7>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<5, dmin<9>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<6, dmin<3>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<6, dmin<5>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<6, dmin<7>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<6, dmin<9>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<7, dmin<3>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<7, dmin<5>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<7, dmin<7>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  cyclic::primitive_bch<7, dmin<9>, normalized_min_sum_tag<50, std::ratio<8, 10>>>(),
  
  cyclic::primitive_bch<5, dmin<3>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<5, dmin<5>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<5, dmin<7>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<5, dmin<9>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<6, dmin<3>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<6, dmin<5>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<6, dmin<7>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<6, dmin<9>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<7, dmin<3>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<7, dmin<5>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<7, dmin<7>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  cyclic::primitive_bch<7, dmin<9>, offset_min_sum_tag<50, std::ratio<1, 100>>>(),
  
  cyclic::primitive_bch<5, dmin<3>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<5>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<9>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<3>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<5>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<7>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<9>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<3>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<5>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<7>, self_correcting_1_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<9>, self_correcting_1_min_sum_tag<50>>(),
  
  cyclic::primitive_bch<5, dmin<3>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<5>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<9>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<3>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<5>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<7>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<9>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<3>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<5>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<7>, self_correcting_2_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<9>, self_correcting_2_min_sum_tag<50>>(),
  
  cyclic::primitive_bch<5, dmin<3>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<5>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<7>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<5, dmin<9>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<3>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<5>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<7>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<6, dmin<9>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<3>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<5>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<7>, normalized_2d_min_sum_tag<50>>(),
  cyclic::primitive_bch<7, dmin<9>, normalized_2d_min_sum_tag<50>>(),
#endif
};

static void usage() {
  int index = 0;
  std::cout << "--algorithm <num>"
            << "  "
            << "Choose algorithm:" << std::endl;

  for (const auto &decoder : decoders)
    std::cout << "  [" << index++ << "] " << decoder.to_string() << std::endl;

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

int main() {
  for (const auto &decoder : decoders) {
    std::cout << decoder.to_string() << std::endl;
  }
  thread_pool p;
  for (const auto &decoder : decoders)
    p.push(awgn_simulation(decoder));
}

