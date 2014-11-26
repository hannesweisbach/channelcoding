#include <vector>
#include <tuple>
#include <iostream>
#include <utility>
#include <iterator>
#include <algorithm>
#include <string>
#include <ctgmath>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <functional>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>
#include <cstdlib>
#include <cctype>

#include "simulation.h"

#include "codes/bch.h"
#include "codes/rs.h"

static std::vector<decoder> decoders{
  cyclic::primitive_bch<5, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<9>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<6, dmin<9>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<3>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<5>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<7, dmin<9>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<3>, min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<5>, min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<9>, min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<3>, min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<5>, min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<7>, min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<9>, min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<3>, min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<5>, min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<7>, min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<9>, min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<3>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<5, dmin<5>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<5, dmin<7>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<5, dmin<9>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<6, dmin<3>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<6, dmin<5>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<6, dmin<7>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<6, dmin<9>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<7, dmin<3>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<7, dmin<5>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<7, dmin<7>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<7, dmin<9>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<5, dmin<3>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<5, dmin<5>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<5, dmin<7>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<5, dmin<9>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<6, dmin<3>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<6, dmin<5>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<6, dmin<7>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<6, dmin<9>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<7, dmin<3>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<7, dmin<5>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<7, dmin<7>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<7, dmin<9>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<5, dmin<3>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<5>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<9>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<3>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<5>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<7>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<9>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<3>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<5>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<7>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<9>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<3>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<5>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<9>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<3>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<5>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<7>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<9>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<3>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<5>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<7>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<9>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<3>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<5>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<9>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<3>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<5>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<7>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<6, dmin<9>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<3>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<5>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<7>, normalized_2d_min_sum_tag<50> >(),
  cyclic::primitive_bch<7, dmin<9>, normalized_2d_min_sum_tag<50> >(),
};

class simulation_factory {
  enum class type {
    AWGN,
    BITFLIP
  };
  enum type type_;

public:
  explicit simulation_factory(const std::string &which = "awgn") {
    if (which == "awgn")
      type_ = type::AWGN;
    else if (which == "bitflip") {
      type_ = type::BITFLIP;
    } else {
      std::ostringstream os;
      os << "Unkown simulation type; " << which;
      throw std::runtime_error(os.str());
    }
  }

  std::function<void(void)> simulation(const class decoder &decoder) const {
    switch (type_) {
    case type::AWGN:
      return awgn_simulation(decoder);
    case type::BITFLIP:
      return bitflip_simulation(decoder);
    }
  }
};

std::string to_lower(const char *s) {
  std::string s_(s);
  std::transform(std::cbegin(s_), std::cend(s_), std::begin(s_), ::tolower);
  return s_;
}

std::string &to_lower(std::string &s) {
  std::transform(std::cbegin(s), std::cend(s), std::begin(s), ::tolower);
  return s;
}

static std::unordered_set<std::string> names;
static std::unordered_set<unsigned> distances;
static std::unordered_set<unsigned> powers;

static std::unordered_multimap<std::string, const decoder &> decoders_by_name;
static std::unordered_multimap<unsigned, const decoder &> decoders_by_distance;
static std::unordered_multimap<unsigned, const decoder &> decoders_by_power;

static void init() {
  for (const auto &decoder : decoders) {
    std::string str = decoder.to_string();
    size_t start;
    size_t end;
    start = str.find_last_of("-");
    std::string name(str.substr(start + 1));
    to_lower(name);
    decoders_by_name.emplace(name, decoder);
    names.emplace(name);

    end = str.find_last_of(")");
    start = str.find_last_of(" ", end);
    auto distance =
        strtoul(str.substr(start + 1, end - start - 1).data(), nullptr, 0);
    decoders_by_distance.emplace(distance, decoder);
    distances.emplace(distance);

    start = str.find_first_of("(");
    end = str.find_first_of(",");
    std::string p(str.substr(start + 1, end - start - 1));
    auto power = std::log2(strtoul(p.data(), nullptr, 0) + 1);

    decoders_by_power.emplace(power, decoder);
    powers.emplace(power);
  }
}

static void usage() {
  std::cout << "--simulation [awgn|bitflip]"
            << "  "
            << "Choose the simulation to run. The default is AWGN."
            << std::endl;
  std::cout << "--algorithm <name>         "
            << "  "
            << "Choose algorithm:" << std::endl;
  for (const auto &algorithm : names)
    std::cout << "  " << algorithm << std::endl;

  std::cout << "--k <num>                  "
            << "  "
            << "Choose code length n = 2^k - 1;" << std::endl;
  for (const auto &k : powers)
    std::cout << "  " << k << std::endl;

  std::cout << "--dmin <num>               "
            << "  "
            << "Choose dmin of the code:" << std::endl;
  for (const auto &dmin : distances)
    std::cout << "  " << dmin << std::endl;

  std::cout << "--seed <num>               "
            << "  "
            << "Set seed of the random number generator. The" << std::endl;
  std::cout << "default is 0." << std::endl;
  std::cout << "--seed-time                "
            << "  "
            << "Use the current time as seed for the random number";
  std::cout << std::endl << "generator." << std::endl;

  std::cout << std::endl;
  std::cout << "algorithm, k, and dmin can be specified multiple times."
            << std::endl;
  std::cout << "For all other options, giving them multiple times results in "
               "the last value" << std::endl << "being used." << std::endl;
  exit(-1);
}

template <typename T> static T convert(const std::string &v) {
  return T(v);
};
template <> std::string convert<std::string>(const std::string &v) {
  return v;
};
template <> unsigned int convert<unsigned int>(const std::string &v) {
  return std::stoul(v);
};

template <typename T>
static void process(const char *option, std::unordered_set<T> &options,
                    const std::unordered_set<T> &reference,
                    const std::string &name) {
  std::string option_string(to_lower(optarg));
  T option_value = convert<T>(option_string);
  if (option_string == "all") {
    options = reference;
  } else if (reference.find(option_value) != std::cend(reference))
    options.emplace(option_value);
  else {
    usage();
    std::ostringstream os;
    os << "Unknown " << name << " '" << option << "'";
    throw std::runtime_error(os.str());
  }
}

static std::tuple<std::unordered_set<std::string>, std::unordered_set<unsigned>,
                  std::unordered_set<unsigned>, simulation_factory, uint64_t>
parse_options(const int argc, char *const argv[]) {

  std::unordered_set<std::string> algorithms;
  std::unordered_set<unsigned> k;
  std::unordered_set<unsigned> dmin;
  simulation_factory factory("awgn");
  uint64_t seed = 0;

  while (1) {
    static struct option options[] = {
      { "simulation", required_argument, nullptr, 'i' },
      { "algorithm", required_argument, nullptr, 'a' },
      { "k", required_argument, nullptr, 'k' },
      { "dmin", required_argument, nullptr, 'd' },
      { "seed", required_argument, nullptr, 's' },
      { "seed-time", no_argument, nullptr, 't' },
      { nullptr, 0, nullptr, 0 },
    };

    int option_index = 0;
    int c = getopt_long_only(argc, argv, "", options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'i': {
      std::string s(to_lower(optarg));
      if (s == "awgn" || s == "bitflip") {
        factory = simulation_factory(s);
      } else {
        std::cerr << "Don't know the simulation type '" << s << "'"
                  << std::endl;
        usage();
      }
    } break;
    case 'a':
      process(optarg, algorithms, names, "algorithm");
      break;
    case 'k':
      process(optarg, k, powers, "k");
      break;
    case 'd':
      process(optarg, dmin, distances, "dmin");
      break;
    case 's':
      seed = strtoull(optarg, nullptr, 0);
      break;
    case 't':
      seed =
          std::chrono::high_resolution_clock::now().time_since_epoch().count();
      break;
    default:
      std::cerr << "Unkown argument: " << c << " " << std::endl;
      usage();
    }
  }

  bool fail = false;

  if (algorithms.empty()) {
    std::cerr << "Algorithm not set." << std::endl;
    fail = true;
  }

  if (k.empty()) {
    std::cerr << "k not set." << std::endl;
    fail = true;
  }

  if (dmin.empty()) {
    std::cerr << "dmin not set." << std::endl;
    fail = true;
  }

  if (fail)
    usage();

  return std::make_tuple(algorithms, k, dmin, factory, seed);
}

int main(int argc, char *const argv[]) {
  std::unordered_set<std::string> algorithms;
  std::unordered_set<unsigned> k;
  std::unordered_set<unsigned> dmin;
  simulation_factory factory;
  uint64_t seed;
  thread_pool p;

  init();
  std::tie(algorithms, k, dmin, factory, seed) = parse_options(argc, argv);

  std::set<const decoder *> chosen_names;
  std::set<const decoder *> chosen_power;
  std::set<const decoder *> chosen_distance;

  for (const auto &algorithm : algorithms) {
    auto iters = decoders_by_name.equal_range(algorithm);
    std::for_each(iters.first, iters.second,
                  [&](const auto &i) { chosen_names.insert(&i.second); });
  }

  for (const auto &power : k) {
    auto iters = decoders_by_power.equal_range(power);
    std::for_each(iters.first, iters.second,
                  [&](const auto &i) { chosen_power.insert(&i.second); });
  }

  for (const auto &distance : dmin) {
    auto iters = decoders_by_distance.equal_range(distance);
    std::for_each(iters.first, iters.second,
                  [&](const auto &i) { chosen_distance.insert(&i.second); });
  }

  std::set<const decoder *> tmp;
  std::set_intersection(std::cbegin(chosen_names), std::cend(chosen_names),
                        std::cbegin(chosen_power), std::cend(chosen_power),
                        std::inserter(tmp, std::end(tmp)));

  std::set<const decoder *> chosen;
  std::set_intersection(
      std::cbegin(tmp), std::cend(tmp), std::cbegin(chosen_distance),
      std::cend(chosen_distance), std::inserter(chosen, std::end(chosen)));

  for (const auto &decoder : chosen) {
    p.push(factory.simulation(*decoder));
  }
}

