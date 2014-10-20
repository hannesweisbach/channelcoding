#pragma once

#include <iostream> 
#include <vector>
#include <algorithm>

#include "gf.h"
#include "center.h"
#include "util.h"

static void protocol_euklid(const std::vector<gf_polynomial> &r,
                            const std::vector<gf_polynomial> &w) {
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

template <typename Field>
gf_polynomial euklid(const Field &field, const unsigned fk,
                     const std::vector<gf_element> &syndromes,
                     const std::vector<gf_element> &erasures =
                         std::vector<gf_element>()) {
  const unsigned dmin = 2 * fk + 1;
  const unsigned max = (2 * fk + erasures.size()) / 2;
  
  gf_polynomial u(field, { 1 });
  std::vector<gf_polynomial> w;
  std::vector<gf_polynomial> r;
  gf_polynomial s(field, syndromes);

  for (const auto &erasure : erasures)
    u *= gf_polynomial(field, { field.one(), erasure });

  r.push_back(s * u);
  r.push_back(gf_polynomial(field));
  std::fill_n(std::back_inserter(r.back()), dmin - 1, field.zero());
  r.back().push_back(field.one());

  w.push_back(u);
  w.push_back(gf_polynomial(field, { 0 }));

  for (int i = 1; r.back().degree() >= max; i++) {
    auto next = r.at(i - 1) % r.at(i);
    auto q = r.at(i - 1) / r.at(i);
    w.push_back(w.at(i - 1) + q * w.at(i));
    r.push_back(next);
  }

  // protocol_euklid(r, w);
  if (w.back().at(0) == field.zero())
    throw decoding_failure("Cannot invert last element");

  return w.back() * w.back().at(0).inverse();
}

template <typename Field>
gf_polynomial berlekamp_massey(const Field &field,
                               const std::vector<gf_element> &syndromes,
                               const std::vector<gf_element> &erasures =
                                   std::vector<gf_element>()) {
#if 0
  std::cout << std::endl << "Calculating Λ(x) with BMA:" << std::endl;
#endif

  const auto fk = syndromes.size() / 2;
  const auto rho = erasures.size();
  gf_polynomial lambda(field, { 1 });
  gf_polynomial b(field);
  int l = erasures.size();

  for (const auto &erasure : erasures)
    lambda *= gf_polynomial(field, { field.one(), erasure });

  b = lambda;

#if 0
  std::cout << "b = " << b << ", lambda: " << lambda << std::endl;
#endif
  for (int i = rho; i < 2 * fk; i++) {
    /* b = b * x; */
    // b *= gf_polynomial(field, { 0, 1 });
    b.push_front(gf_element(field, 0));
    /* d = si + \sigma_j=1^l lambda_j * s_i-j; */
    const auto delta =
        std::inner_product(std::cbegin(lambda) + 1, std::cbegin(lambda) + l + 1,
                           std::crend(syndromes) - i, syndromes.at(i));

    if (delta) {
      auto t = lambda + b * delta;
//      std::cout << "T(x) = " << t << std::endl;
      if (2 * l <= i + rho) {
        b = lambda * delta.inverse();
        l = i + rho - l + 1;
      }
      lambda = t;
    }
#if 0
    std::cout << i << " | " << delta << " | " << lambda << " | " << l << " | "
              << b << " " << std::endl;
#endif
  }

  //std::cout << "lambda: " << lambda << std::endl;
  return lambda;
}

