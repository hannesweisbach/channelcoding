#include <iostream>
#include <cstdlib>
#include <string>

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
  unsigned l = 0;
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
      l = static_cast<unsigned>(std::stoul(optarg));
      break;
    case 's':
      seed = std::stoull(optarg);
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
