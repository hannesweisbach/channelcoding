#include <algorithm>

#include "bch.h"
#include "iterative.h"

unsigned test_hard(const unsigned errors, const matrix<int>& H) {
  std::vector<std::vector<int>> wrong_reconstruction;
  std::vector<std::vector<int>> reconstruction_failure;
  const size_t length = 31;
  std::vector<int> b;
  
  std::fill_n(std::back_inserter(b), length - errors, 0);
  std::fill_n(std::back_inserter(b), errors, 1);

  do {
    std::vector<float> x;
    x.reserve(length);

    std::transform(std::cbegin(b), std::cend(b), std::back_inserter(x),
                   [](const auto &bit) { return -2 * bit + 1; });
    
    try {
      auto result = min_sum<1000, float>(H, x);
      const auto& b_corr = std::get<0>(result);
      auto it = std::find_if(std::cbegin(b_corr), std::cend(b_corr),
                             [](const auto &bit) { return bit != 0; });
      if (it != std::cend(b_corr)) {
        wrong_reconstruction.push_back(b);
#if 1
        std::cout << "b  = ";
        for (const auto &bit : b)
          std::cout << bit;
        std::cout << " got wrongly corrected to:" << std::endl;
        std::cout << "b_ = ";
        for (const auto &bit : b_corr)
          std::cout << bit;
        std::cout << " after " << std::get<2>(result) << " iterations"
                  << std::endl;
#endif
      }
    }
    catch (...) {
      reconstruction_failure.push_back(b);
#if 0
      std::cout << "Reconstruction failure for b = ";
      for (const auto &bit : b)
        std::cout << bit;
      std::cout << std::endl;
#endif
    }
  } while(std::next_permutation(std::begin(b), std::end(b)));

  std::cout << "Experienced " << wrong_reconstruction.size()
            << " wrong reconstructions and " << reconstruction_failure.size()
            << " reconstruction failures for " << errors << " errors"
            << std::endl;

  return reconstruction_failure.size();
}

int main() {
  const unsigned fk = 3;

  bch code(5, 0x25, 7);

  for (int err = 0; err < fk + 1; err++) {
    test_hard(err, code.H());
  }
}

