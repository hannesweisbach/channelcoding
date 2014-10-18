#include <algorithm>
#include <fstream>
#include <chrono>

#include "bch.h"
#include "iterative.h"
#include "util.h"

using decoder_t =
    std::function<std::tuple<std::vector<int>, std::vector<float>, unsigned>(
        std::vector<float>)>;

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
    catch (...) {
      failures++;
    }

  } while (std::next_permutation(std::begin(b), std::end(b)));

  return (double)failures / error_patterns;
}

std::tuple<std::vector<int>, std::vector<float>, unsigned>
pzg_wrapper(const bch &code, const std::vector<float> &b) {
  std::vector<int> hard;
  hard.reserve(b.size());
  std::transform(std::cbegin(b), std::cend(b), std::back_inserter(hard),
                 [](const auto &bit) { return bit < 0; });
  auto b_corr = code.correct_bm(hard);
  return std::make_tuple(b_corr, std::vector<float>(), 0);
}

int main() {
  constexpr size_t max_iterations = 50;
  const unsigned max_fk = 7;

  bch code(5, 0x25, 7);

  std::vector<std::pair<std::string, decoder_t> > algorithms;
  algorithms.emplace_back(
      "pzg", std::bind(pzg_wrapper, std::cref(code), std::placeholders::_1));
  algorithms.emplace_back("ms", std::bind(min_sum<max_iterations, float, float>,
                                          std::cref(code.H()),
                                          std::placeholders::_1));
  algorithms.emplace_back("nms", std::bind(nms<max_iterations, float, float>,
                                           std::cref(code.H()),
                                           std::placeholders::_1, 0.9f));
  algorithms.emplace_back(
      "scms1", std::bind(scms1<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));
  algorithms.emplace_back(
      "scms2", std::bind(scms2<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));

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

