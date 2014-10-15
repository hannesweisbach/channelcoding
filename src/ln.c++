#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>

#include "util.h"
#include "bch.h"
#include "iterative.h"

static float sigma(const float eb_n0, const float R) {
  return sqrt(1.0f / (2 * R * pow(10, eb_n0 / 10.0)));
}

void generate_hard(const size_t samples, const float eb_n0_max,
                   const float eb_n0_step = 0.1f) {

  bch code(5, 0x25, 7);
  constexpr int fk = 3;
  constexpr size_t length = 31;
  constexpr float R = 16.0/length;
  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  std::vector<float> b(length);

  for (float eb_n0 = 0; eb_n0 < eb_n0_max; eb_n0 += eb_n0_step) {
    std::normal_distribution<float> random(1.0, sigma(eb_n0, R));
    auto normal_gen = std::bind(std::ref(random), std::ref(generator));
    
    unsigned failures = 0;
    unsigned fk_corr = 0;
    unsigned bit_errors = 0;

    for (size_t sample = 0; sample < samples; sample++) {
      std::generate(std::begin(b), std::end(b), normal_gen);
      try {
        auto b_corr = code.correct_peterson(hard_decision(b));
        auto wrong_bits =
            std::count_if(std::cbegin(b_corr), std::cend(b_corr),
                          [](const auto &bit) { return bit != 0; });
        if (wrong_bits) {
          failures++;
          bit_errors += wrong_bits;
        } else {
          auto d = std::count_if(std::cbegin(b), std::cend(b),
                                 [](const auto &bit) { return bit < 0; });
          if (d > fk)
            fk_corr++;
        }
      }
      catch (...) {
#if 0
        for (auto e : b)
          std::cout << e << ", ";
        std::cout << std::endl;
#endif
        failures++;
      }
    }
    /* word error rate */
    std::cout << std::scientific;
    std::cout << std::setprecision(4) << eb_n0 << " ";
    std::cout << std::setprecision(9) << failures / (float)samples << " ";
    std::cout << std::setprecision(9) << fk_corr / (float)samples << " ";
    std::cout << std::setprecision(9) << bit_errors / (float)(samples * length);
    std::cout << std::endl;
  }

}

void generate_diagram(const matrix<int> &H, const size_t samples,
                      const float eb_n0_max, const float eb_n0_step = 0.1f) {
  constexpr int fk = 3;
  constexpr size_t length = 31;
  constexpr float R = 16.0/length;
  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  
  std::vector<float> b(length);

  for (float eb_n0 = 0; eb_n0 < eb_n0_max; eb_n0 += eb_n0_step) {
    std::normal_distribution<float> random(1.0, sigma(eb_n0, R));
    auto normal_gen = std::bind(std::ref(random), std::ref(generator));
    
    unsigned failures = 0;
    unsigned fk_corr = 0;
    unsigned bit_errors = 0;

    for (size_t sample = 0; sample < samples; sample++) {
      std::generate(std::begin(b), std::end(b), normal_gen);
      try {
        auto result = min_sum<200, float>(H, b);
        const auto& b_corr = std::get<0>(result);
        auto wrong_bits =
            std::count_if(std::cbegin(b_corr), std::cend(b_corr),
                          [](const auto &bit) { return bit != 0; });
#if 0
        std::cout << wrong_bits << ": ";  
        for (auto e : b_corr)
            std::cout << e << ", ";
          std::cout << std::endl;
#endif
        if (wrong_bits) {
          failures++;
          bit_errors += wrong_bits;
        } else {
          auto d = std::count_if(std::cbegin(b), std::cend(b),
                                 [](const auto &bit) { return bit < 0; });
          if (d > fk)
            fk_corr++;
        }
      }
      catch (...) {
#if 0
        for (auto e : b)
          std::cout << e << ", ";
        std::cout << std::endl;
#endif
        failures++;
      }
    }
    /* word error rate */
    std::cout << std::scientific;
    std::cout << std::setprecision(4) << eb_n0 << " ";
    std::cout << std::setprecision(9) << failures / (float)samples << " ";
    std::cout << std::setprecision(9) << fk_corr / (float)samples << " ";
    std::cout << std::setprecision(9) << bit_errors / (float)(samples * length);
    std::cout << std::endl;
  }
}


int main() {
  const unsigned errors = 2;
  const unsigned N = 10000000;
  std::mt19937 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<int> distribution(0,30);

  bch code(5, 0x25, 7);

  std::cout << code.H() << std::endl;
  std::vector<int> a(31, 0);

  unsigned failures = 0;

  generate_hard(10000, 10);
  return 0;
  
  generate_diagram(code.H(), 10000, 10);

  return 0;
#if 0
  for (int err = 0; err < 4; err++) {
    test_hard(err, code.H());
  }
#endif

  std::vector<int> test({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 });
  std::vector<int> test2({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, });
  std::vector<float> x;
  for (const auto &e : test2) {
    x.push_back(e ? -1 : 1);
    std::cout << x.back() << " ";
  }
  std::cout << std::endl;

  min_sum<1000, float>(code.H(), x);

  return 0;

  for (size_t i = 0; i < N; i++) {
    auto b = a;

    for (unsigned error = 0; error < errors; error++)
      b.at(distribution(generator)) ^= 1;

#if 0
    std::cout << "code word: ";
    for (const auto &e : a)
      std::cout << e << " ";
    std::cout << std::endl;

    std::cout << "recv word: ";
    for (const auto &e : b)
      std::cout << e << " ";
    std::cout << std::endl;
#endif

    std::vector<float> x;
    for (const auto &e : b)
      x.push_back(e ? -1 : 1);
#if 0
    for (auto e : x)
      std::cout << e << " ";
    std::cout << std::endl;
#endif
    
    try {
      auto result = min_sum<10000, float>(code.H(), x);
    }
    catch (...) {
      failures++;
      std::cout << failures << " in " << i << " tries" << std::endl;
    }
  }
}
