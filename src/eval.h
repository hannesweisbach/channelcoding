#pragma once

#include <ostream>
#include <random>
#include <string>

#include "bch.h"
#include "iterative.h"

using decoder_t =
    std::function<std::tuple<std::vector<int>, std::vector<float>, unsigned>(
        std::vector<float>)>;

uint64_t primitive_polynomial(const unsigned degree);

std::tuple<std::vector<int>, int, int, uint64_t>
parse_options(const int argc, char *const argv[]);

inline size_t iterations(const double ber) {
  return std::min(1000 * 1 / ber, 10e6);
  // return std::min(base_trials * pow(10, eb_no / 2), 10e6);
};

std::tuple<std::vector<int>, std::vector<float>, unsigned>
pzg_wrapper(const bch &code, const std::vector<float> &b);

std::tuple<std::vector<int>, std::vector<float>, unsigned>
bm_wrapper(const bch &code, const std::vector<float> &b);

std::tuple<std::vector<int>, std::vector<float>, unsigned>
uncoded(const std::vector<float> &b);

template <unsigned max_iterations>
std::vector<std::pair<std::string, decoder_t> >
get_algorithms(const bch &code) {
  std::vector<std::pair<std::string, decoder_t> > algorithms;
  algorithms.emplace_back("ms", std::bind(min_sum<max_iterations, float, float>,
                                          std::cref(code.H()),
                                          std::placeholders::_1));
  algorithms.emplace_back(
      "scms1", std::bind(scms1<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));
  algorithms.emplace_back(
      "scms2", std::bind(scms2<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1));

  algorithms.emplace_back(
      "pzg", std::bind(pzg_wrapper, std::cref(code), std::placeholders::_1));
  algorithms.emplace_back(
      "bm", std::bind(bm_wrapper, std::cref(code), std::placeholders::_1));
  algorithms.emplace_back("nms", std::bind(nms<max_iterations, float, float>,
                                           std::cref(code.H()),
                                           std::placeholders::_1, 0.915f));
  algorithms.emplace_back("oms", std::bind(oms<max_iterations, float, float>,
                                           std::cref(code.H()),
                                           std::placeholders::_1, 0.029f));
  algorithms.emplace_back("2d-nms",
                          std::bind(nms_2d<max_iterations, float, float>,
                                    std::cref(code.H()), std::placeholders::_1,
                                    0.97f, 0.92f));
  return algorithms;
}

class eval_object {
public:
  const float eb_no;
  const unsigned reconstruction_failures = 0;
  const unsigned reconstruction_errors = 0;
  const unsigned fk_corr = 0;
  const unsigned bit_errors = 0;
  const size_t samples = 0;
  const unsigned n = 0;

  double wer() const;
  double ber() const;

  static std::string header();
  friend std::ostream &operator<<(std::ostream &os, const eval_object &);
};

eval_object evaluate(std::mt19937_64 &generator, const decoder_t &decoder,
                     const size_t samples, const float eb_n0, const unsigned n,
                     const unsigned l, const unsigned fk);

inline auto build_simulator(std::mt19937_64 &generator, const bch &code) {
  auto parameters = code.parameters();
  return std::bind(evaluate, generator, std::placeholders::_1,
                   std::placeholders::_2, std::placeholders::_3,
                   std::get<0>(parameters), std::get<1>(parameters),
                   std::get<2>(parameters));
}
