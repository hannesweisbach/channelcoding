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

  const float alpha_start = (argc > 1) ? strtof(argv[1], nullptr) : 0.1f;
  const float alpha_max = (argc > 2) ? strtof(argv[2], nullptr) : 1.0f;
  const float alpha_step = (argc > 3) ? strtof(argv[3], nullptr) : 0.01f;

  const float beta_start = (argc > 4) ? strtof(argv[4], nullptr) : 0.1f;
  const float beta_max = (argc > 5) ? strtof(argv[5], nullptr) : 1.0f;
  const float beta_step = (argc > 6) ? strtof(argv[6], nullptr) : 0.01f;

  const float eb_no = (argc > 7) ? strtof(argv[7], nullptr) : 0.0f;

  std::cout << "Calculating for eb_no: " << eb_no << std::endl;
  std::cout << "alpha: " << alpha_start << " to " << alpha_max << " in "
            << alpha_step << " steps" << std::endl;
  std::cout << "beta:  " << beta_start << " to " << beta_max << " in "
            << beta_step << " steps" << std::endl;

  std::mt19937_64 generator;
  const auto seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

  bch code(5, 0x25, 7);
  // bch code(6, 0x45, 7);

  auto parameters = code.parameters();
  auto simulator = std::bind(evaluate, generator, std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             std::get<0>(parameters), std::get<1>(parameters),
                             std::get<2>(parameters));

  std::ostringstream os;
  os << "logs/ebno_" << eb_no;
  std::string fname = os.str();
  if (file_exists(fname)) {
    std::cout << "File " << fname << " exists. Aborting.";
    return -1;
  }

  std::ofstream file(fname, std::ofstream::out);
  std::cout << "Writing to file " << fname << std::endl;
  file << "#eb_no " << eb_no << std::endl;
  file << "#alpha: " << alpha_start << " " << alpha_max << " " << alpha_step
       << " in rows" << std::endl;
  file << "#beta : " << beta_start << " " << beta_max << " " << beta_step
       << " in columns" << std::endl;

  file << std::scientific;

  double wer = 0.5;
  for (float alpha = alpha_start; alpha < alpha_max; alpha += alpha_step) {
    for (float beta = beta_start; beta < beta_max; beta += beta_step) {
      generator.seed(seed);

      const size_t N = iterations(wer);

      std::cout << "alpha = " << alpha << " beta = " << beta << " N = " << N;
      std::cout.flush();

      auto start = std::chrono::high_resolution_clock::now();
      auto f =
          std::bind(nms_2d<max_iterations, float, float>, std::ref(code.H()),
                    std::placeholders::_1, alpha, beta);

      wer = simulator(f, N, eb_no).wer();
      file << wer << " ";

      auto end = std::chrono::high_resolution_clock::now();
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

      std::cout << seconds << " s" << std::endl;
    }

    file << std::endl;
  }
}
