#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <cstdlib>

#include "eval.h"
#include "util.h"

int main(int argc, const char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;

  float beta_start = (argc > 1) ? strtof(argv[1], nullptr) : 0.1f;
  float beta_max = (argc > 2) ? strtof(argv[2], nullptr) : 1.0f;
  float beta_step = (argc > 3) ? strtof(argv[3], nullptr) : 0.01f;

  std::cout << "Calculating from " << beta_start << " to " << beta_max << " in "
            << beta_step << " steps" << std::endl;

  std::mt19937_64 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  int k = 5;
  bch code(k, primitive_polynomial(k), 7);
  auto simulator = build_simulator(generator, code);

  const char *fname;
  std::ostringstream os;
  auto parameters = code.parameters();
  os << "logs/beta_" << std::setprecision(3) << std::fixed << beta_start << "_"
     << std::get<0>(parameters) << "_" << std::get<1>(parameters) << "_"
     << std::get<2>(parameters) << ".dat";

  if (file_exists(os.str())) {
    std::cerr << "File " << fname << " already exists." << std::endl;
    return -1;
  }
  fname = os.str().c_str();

  std::ofstream file(fname, std::ofstream::out);
  std::cout << "Writing to file " << fname << std::endl;

  /* header */
  file << "eb_no ";
  for (float beta = beta_start; beta < beta_max; beta += beta_step) {
    std::ostringstream col_name;
    col_name << "ber_" << std::fixed << std::setprecision(3) << beta << " ";
    file << col_name.str();
  }
  file << std::endl;
  file << std::scientific;

  for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
    file << eb_no << " ";

    std::cout << "Calculating for E_b/N_0 = " << eb_no << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (float beta = beta_start; beta < beta_max; beta += beta_step) {
      auto f = std::bind(oms<max_iterations, float, float>, std::ref(code.H()),
                         std::placeholders::_1, beta);
      const size_t N = iterations(eb_no);
      file << simulator(f, N, eb_no).ber() << " ";
      std::cout << beta << " " << N << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << seconds << " s per iteration" << std::endl;
    file << std::endl;
  }
}
