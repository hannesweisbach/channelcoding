#include <random>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <cstdlib>

#include "util.h"
#include "eval.h"

int main(int argc, const char *const argv[]) {
  constexpr size_t max_iterations = 50;
  constexpr float eb_no_max = 10;
  constexpr float eb_no_step = 0.1f;

  float alpha_start = (argc > 1) ? strtof(argv[1], nullptr) : 0.1f;
  float alpha_max = (argc > 2) ? strtof(argv[2], nullptr) : 1.0f;
  float alpha_step = (argc > 3) ? strtof(argv[3], nullptr) : 0.01f;

  std::cout << "Calculating from " << alpha_start << " to " << alpha_max
            << " in " << alpha_step << " steps" << std::endl;

  std::mt19937_64 generator(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  bch code(5, 0x25, 7);
  // bch code(6, 0x45, 7);

  auto parameters = code.parameters();
  auto simulator = std::bind(evaluate, generator, std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             std::get<0>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));

  std::ostringstream os;
  os << "alpha_" << alpha_start << "_" << std::get<0>(parameters) << "_"
     << std::get<1>(parameters) << "_" << std::get<2>(parameters) << ".dat";
  std::string fname(os.str());
  if (file_exists(fname)) {
    std::cerr << "File " << fname << " already exists." << std::endl;
    return -1;
  }

  std::ofstream file(fname, std::ofstream::out);
  std::cout << "Writing to file " << fname << std::endl;

  /* header */
  file << "eb_no ";
  for (float alpha = alpha_start; alpha < alpha_max; alpha += alpha_step) {
    std::ostringstream col_name;
    col_name << "ber_" << alpha * 1000 << " ";
    file << col_name.str();
  }
  file << std::endl;
  file << std::scientific;

  for (double eb_no = 0; eb_no < eb_no_max; eb_no += eb_no_step) {
    file << eb_no << " ";

    std::cout << "Calculating for E_b/N_0 = " << eb_no << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (float alpha = alpha_start; alpha < alpha_max; alpha += alpha_step) {
#if 0
      auto f = std::bind(nms<max_iterations, float, float>, std::ref(code.H()),
                         std::placeholders::_1, alpha);
#else
      auto f = std::bind(scms1_nms<max_iterations, float, float>,
                         std::cref(code.H()), std::placeholders::_1, alpha);
#endif
      const size_t N = iterations(eb_no);
      file << simulator(f, N, eb_no).ber() << " ";
      std::cout << eb_no << " " << alpha << " " << N << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    std::cout << seconds << " s per iteration" << std::endl;
    file << std::endl;
  }
}
