#include <iostream>
#include <cstdlib>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>

#include "codes/uncoded.h"
#include "simulation/simulation.h"

[[noreturn]] static void usage() {
  std::cout << "--l <num>        "
            << "  "
            << "Choose code length l" << std::endl;
  std::cout << "--seed <num>     "
            << "  "
            << "Set seed of the random number generator." << std::endl;
  std::exit(EXIT_FAILURE);
}

int main(int argc, char *const argv[]) {
  unsigned long l = 0;
  uint64_t seed = 0;

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

  awgn_simulation(uncoded(l), 0.5, seed)();
}
