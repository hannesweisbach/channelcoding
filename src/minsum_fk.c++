#include <algorithm>
#include <fstream>
#include <chrono>

#include "eval.h"

double test_hard(const unsigned errors, const decoder_t &decoder) {
  const unsigned length = 31;
  size_t error_patterns = 0;
  size_t failures = 0;
  std::vector<int> b;

  std::fill_n(std::back_inserter(b), length - errors, 0);
  std::fill_n(std::back_inserter(b), errors, 1);

  do {
    std::vector<float> x;
    x.reserve(length);

    std::transform(std::cbegin(b), std::cend(b), std::back_inserter(x),
                   [](const auto &bit) { return -2 * bit + 1; });

    error_patterns++;
    try {
      auto result = decoder(x);
      const auto &b_corr = std::get<0>(result);
      auto it = std::find_if(std::cbegin(b_corr), std::cend(b_corr),
                             [](const auto &bit) { return bit != 0; });
      if (it != std::cend(b_corr)) {
        failures++;
      }
    }
    catch (const decoding_failure &) {
      failures++;
    }

  } while (std::next_permutation(std::begin(b), std::end(b)));

  return (double)failures / error_patterns;
}

int main() {
  constexpr size_t max_iterations = 50;
  const unsigned max_fk = 7;

  //bch code(5, 0x25, 7);
  bch code(6, 0x43, 7);

  auto algorithms = get_algorithms<max_iterations>(code);

  char fname[] = "fk.XXXXXX";
  mktemp(fname);
  std::ofstream file(fname, std::ofstream::out);

  std::cout << "Writing to file " << fname << std::endl;

  file << "err ";
  for (const auto &algorithm : algorithms) {
    file << algorithm.first << " ";
  }
  file << std::endl;

  for (int err = 0; err < max_fk; err++) {
    file << err << " ";
    for (const auto &algorithm : algorithms) {
      const auto &name = algorithm.first;
      const auto &decoder = algorithm.second;

      std::cout << "Running " << name << " with " << err << " errors … ";
      std::cout.flush();
      auto start = std::chrono::high_resolution_clock::now();

      file << test_hard(err, decoder) << " ";

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

      std::cout << seconds << " s" << std::endl;
    }

    file << std::endl;
  }
}

