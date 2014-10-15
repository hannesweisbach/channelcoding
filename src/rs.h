#pragma once

#include <vector>

#include "gf.h"

class rs {
public:
  rs(unsigned q, uint64_t modular_polynom, unsigned dmin);

  gf_polynomial correct_pzg(const gf_polynomial &b) const;
  gf_polynomial correct_erasures(const gf_polynomial &b,
                                 const std::vector<gf_element> &erasures) const;
  gf_polynomial correct_bm(const gf_polynomial &b,
                           const std::vector<gf_element> &erasures =
                               std::vector<gf_element>()) const;

private:
  gf_polynomial pzg(const std::vector<gf_element> &syndromes) const;
  gf_polynomial error_values_naive(const std::vector<gf_element> &syndromes,
                                   const std::vector<gf_element> &zeroes) const;
  gf field;

  unsigned n;
  unsigned l;
  unsigned k;
  unsigned dmin;
  unsigned mu;
  unsigned fk;

  gf_polynomial generator;
};
