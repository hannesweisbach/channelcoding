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

std::ostream &operator<<(std::ostream &os, const eval_object &e) {
  os << std::scientific;
  const double samples_ = e.samples;
  os << std::setprecision(12) << e.reconstruction_failures / samples_ << " ";
  os << std::setprecision(12) << e.reconstruction_errors / samples_ << " ";
  const auto wrong_words = e.reconstruction_failures + e.reconstruction_errors;
  os << std::setprecision(12) << wrong_words / samples_ << " ";
  os << std::setprecision(12) << e.fk_corr / samples_ << " ";
  os << std::setprecision(12) << e.bit_errors / (samples_ * e.l);
  return os;
}

double eval_object::ber() const { return bit_errors / (double)(samples * l); }

