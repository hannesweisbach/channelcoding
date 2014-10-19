#include "eval.h"

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

