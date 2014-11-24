#include <vector>

#include "simulation.h"

#include "codes/bch.h"
#include "codes/rs.h"

std::vector<decoder> decoders{
  cyclic::primitive_bch<5, dmin<7>, cyclic::berlekamp_massey_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::peterson_gorenstein_zierler_tag>(),
  cyclic::primitive_bch<5, dmin<7>, cyclic::euklid_tag>(),
  cyclic::primitive_bch<5, dmin<7>, min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>,
                        normalized_min_sum_tag<50, std::ratio<8, 10> > >(),
  cyclic::primitive_bch<5, dmin<7>,
                        offset_min_sum_tag<50, std::ratio<1, 100> > >(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_1_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, self_correcting_2_min_sum_tag<50> >(),
  cyclic::primitive_bch<5, dmin<7>, normalized_2d_min_sum_tag<50> >(),
};

int main() {
  for (const auto &decoder : decoders) {
    std::cout << decoder.to_string() << std::endl;
  }
  thread_pool p;
  for (const auto &decoder : decoders)
    p.push(bitflip_simulation(decoder, 6));
}

