#pragma once

#include "bch.h"
#include "iterative.h"

using decoder_t =
    std::function<std::tuple<std::vector<int>, std::vector<float>, unsigned>(
        std::vector<float>)>;

std::tuple<std::vector<int>, std::vector<float>, unsigned>
pzg_wrapper(const bch &code, const std::vector<float> &b);

std::tuple<std::vector<int>, std::vector<float>, unsigned>
bm_wrapper(const bch &code, const std::vector<float> &b);

template <unsigned max_iterations>
std::vector<std::pair<std::string, decoder_t> >
get_algorithms(const bch &code) {
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
  return algorithms;
}
