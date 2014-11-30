#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <ostream>
#include <algorithm>
#include <iterator>
#include <clocale>
#include <type_traits>

#include "center.h"

template <typename Polynomial>
static void protocol_euklid(const std::vector<Polynomial> &,
                            const std::vector<Polynomial> &, std::false_type) {}

template <typename Polynomial>
static void protocol_euklid(const std::vector<Polynomial> &r,
                            const std::vector<Polynomial> &w, std::true_type) {
  std::setlocale(LC_ALL, "en_US.UTF-8");

  std::vector<std::string> r_;
  std::vector<std::string> w_;
  std::vector<std::string> q_;

  auto poly2str = [](const auto &poly) {
    std::ostringstream os;
    os << poly;
    return os.str();
  };

  std::transform(std::cbegin(r), std::cend(r), std::back_inserter(r_),
                 poly2str);
  std::transform(std::cbegin(w), std::cend(w), std::back_inserter(w_),
                 poly2str);
  std::transform(std::cbegin(r) + 1, std::cend(r), std::cbegin(r),
                 std::back_inserter(q_),
                 [&](const auto &first1, const auto &first2) {
    auto q = first2 / first1;
    return poly2str(q);
  });

  auto mb_strcmp = [](const std::string &lhs, const std::string &rhs) {
    return mbstowcs(nullptr, lhs.c_str(), 0) <
           mbstowcs(nullptr, rhs.c_str(), 0);
  };
  const auto max_w =
      std::max_element(std::cbegin(w_), std::cend(w_), mb_strcmp);
  const auto max_r =
      std::max_element(std::cbegin(r_), std::cend(r_), mb_strcmp);
  const auto max_q =
      std::max_element(std::cbegin(q_), std::cend(q_), mb_strcmp);

  const auto w_len = mbstowcs(nullptr, max_w->c_str(), 0);
  const auto r_len = mbstowcs(nullptr, max_r->c_str(), 0);
  const auto q_len = mbstowcs(nullptr, max_q->c_str(), 0);
  const auto i_len = std::max(r.size() / 3, 2UL);

  std::cout << std::endl << "Protocol EUKLID algorithm:" << std::endl;

  std::cout << std::setw(i_len) << std::right << "i"
            << " | ";
  std::cout << std::setw(r_len) << centered("r(x)") << " | ";
  std::cout << std::setw(q_len) << centered("q(x)") << " | ";
  std::cout << std::setw(w_len) << centered("w(x)") << " |" << std::endl;

  std::cout << std::setfill('-') << std::setw(i_len) << ""
            << "-+-";
  std::cout << std::setw(r_len) << ""
            << "-+-";
  std::cout << std::setw(q_len) << ""
            << "-+-";
  std::cout << std::setw(w_len) << ""
            << "-+";
  std::cout << std::endl << std::setfill(' ');

  for (int i = 1; i < r.size(); i++) {
    std::cout << std::setw(i_len) << std::right << i - 2 << " | ";
    std::cout << std::setw(r_len) << centered(r_.at(i)) << " | ";
    if (i < 2) {
      std::cout << std::setw(q_len) << ""
                << " | ";
    } else {
      std::cout << std::setw(q_len) << centered(q_.at(i - 1)) << " | ";
    }
    std::cout << std::setw(w_len) << centered(w_.at(i)) << " |" << std::endl;
  }

  std::cout << std::setfill('-') << std::setw(i_len) << ""
            << "-+-";
  std::cout << std::setw(r_len) << ""
            << "-+-";
  std::cout << std::setw(q_len) << ""
            << "-+-";
  std::cout << std::setw(w_len) << ""
            << "-+";
  std::cout << std::endl << std::setfill(' ');

  std::cout << "Λ(x)_euklid = " << w.back() << " * " << w.back().at(0)
            << "^-1 = " << w.back() * w.back().at(0).inverse() << std::endl;
  std::cout << "v = " << w.back().degree() - w.at(0).degree() << std::endl
            << std::endl;
}

template <typename Polynomial,
          typename Element = typename Polynomial::coefficient_type>
void protocol_bm(const size_t &, const Element &, const Polynomial &,
                 const size_t &, const Polynomial &, std::false_type) {}

template <typename Polynomial,
          typename Element = typename Polynomial::coefficient_type>
void protocol_bm(const size_t &i, const Element &delta,
                 const Polynomial &lambda, const size_t &l, const Polynomial &b,
                 std::true_type) {
  std::cout << i << " | " << delta << " | " << lambda << " | " << l << " | "
            << b << " " << std::endl;
}
