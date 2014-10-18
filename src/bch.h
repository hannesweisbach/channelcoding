#pragma once

#include <vector>

#include "matrix.h"
#include "gf.h"

std::vector<std::vector<unsigned int> > bch_cycles(unsigned int n);
std::vector<gf_element> chien(const gf &field, gf_polynomial polynomial);

class bch {
public:
  enum class type {
    shortened,
    extended,
    normal
  };

  enum class coding {
    multiplication,
    division
  };

  enum class locator_method {
    PGZ,
    EUKLID,
    BMA
  };
  enum class zeroes {
    brute_force,
    chien
  };
  bch(unsigned q, uint64_t modular_polynomial, unsigned d_e);
  bch(unsigned q, uint64_t modular_polynomial, unsigned d_e, unsigned n);

  std::vector<int> encode_div(const std::vector<int> &b) const;
  
  std::vector<gf_element> syndromes(const gf_polynomial &a) const;
  gf_polynomial pzg(std::vector<gf_element> syndromes) const;

  std::vector<int> correct_peterson(const std::vector<int> &b) const;
  std::vector<int> decode(const std::vector<int> &b) const;
  //decode_mult
  //decode_div
  //encode_mult
  gf_polynomial correct_bm(const gf_polynomial &b,
                           const std::vector<gf_element> &erasures =
                               std::vector<gf_element>()) const;
  std::vector<int> correct_bm(const std::vector<int> &b,
                              const std::vector<gf_element> &erasures =
                                  std::vector<gf_element>()) const;
  const matrix<int> &H() const;

  gf field;
private:

  unsigned n;
  unsigned l;
  unsigned k;
  unsigned dmin;
  unsigned mu;
  unsigned fk;

  enum type type = type::normal;

  gf_polynomial generator;
  gf_polynomial f;
  gf_polynomial h;
  matrix<int> control_matrix;
};

