#include <cmath>
#include <array>
#include <chrono>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <sys/stat.h>

#include "simulation.h"
#include "codes/codes.h"

namespace detail {
inline bool file_exists(const std::string &fname) {
  struct stat buf;
  return (stat(fname.c_str(), &buf) != -1);
}
}

static std::array<double, 131> rates{
  { 0.01,  0.02,  0.03,  0.04,  0.05,  0.06,  0.07,  0.08,  0.09,  0.10,  0.11,
    0.12,  0.13,  0.14,  0.15,  0.16,  0.17,  0.18,  0.19,  0.20,  0.21,  0.22,
    0.23,  0.24,  0.25,  0.26,  0.27,  0.28,  0.29,  0.30,  0.31,  0.32,  0.33,
    0.34,  0.35,  0.36,  0.37,  0.38,  0.39,  0.40,  0.41,  0.42,  0.43,  0.44,
    0.45,  0.46,  0.47,  0.48,  0.49,  0.50,  0.51,  0.52,  0.53,  0.54,  0.55,
    0.56,  0.57,  0.58,  0.59,  0.60,  0.61,  0.62,  0.63,  0.64,  0.65,  0.66,
    0.67,  0.68,  0.69,  0.70,  0.71,  0.72,  0.73,  0.74,  0.75,  0.76,  0.77,
    0.78,  0.79,  0.800, 0.807, 0.817, 0.827, 0.837, 0.846, 0.855, 0.864, 0.872,
    0.880, 0.887, 0.894, 0.900, 0.907, 0.913, 0.918, 0.924, 0.929, 0.934, 0.938,
    0.943, 0.947, 0.951, 0.954, 0.958, 0.961, 0.964, 0.967, 0.970, 0.972, 0.974,
    0.976, 0.978, 0.980, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989,
    0.990, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999 }
};

static std::array<double, 131> limits = {
  { -1.548, -1.531, -1.500, -1.470, -1.440, -1.409, -1.378, -1.347, -1.316,
    -1.285, -1.254, -1.222, -1.190, -1.158, -1.126, -1.094, -1.061, -1.028,
    -0.995, -0.963, -0.928, -0.896, -0.861, -0.827, -0.793, -0.757, -0.724,
    -0.687, -0.651, -0.616, -0.579, -0.544, -0.507, -0.469, -0.432, -0.394,
    -0.355, -0.314, -0.276, -0.236, -0.198, -0.156, -0.118, -0.074, -0.032,
    0.010,  0.055,  0.097,  0.144,  0.188,  0.233,  0.279,  0.326,  0.374,
    0.424,  0.474,  0.526,  0.574,  0.628,  0.682,  0.734,  0.791,  0.844,
    0.904,  0.960,  1.021,  1.084,  1.143,  1.208,  1.275,  1.343,  1.412,
    1.483,  1.554,  1.628,  1.708,  1.784,  1.867,  1.952,  2.045,  2.108,
    2.204,  2.302,  2.402,  2.503,  2.600,  2.712,  2.812,  2.913,  3.009,
    3.114,  3.205,  3.312,  3.414,  3.500,  3.612,  3.709,  3.815,  3.906,
    4.014,  4.115,  4.218,  4.304,  4.425,  4.521,  4.618,  4.725,  4.841,
    4.922,  5.004,  5.104,  5.196,  5.307,  5.418,  5.484,  5.549,  5.615,
    5.681,  5.756,  5.842,  5.927,  6.023,  6.119,  6.234,  6.360,  6.495,
    6.651,  6.837,  7.072,  7.378,  7.864 }
};

decoder::decoder_concept::~decoder_concept() = default;

static constexpr double ebno(const double rate) {
  if (rate <= 0.800) {
    const size_t index = static_cast<size_t>(rate * 100);
    return limits.at(index);
  } else if (rate >= 0.999) {
    return limits.back();
  } else {
    const size_t offset = 80;
    const size_t index = static_cast<size_t>(std::distance(
        std::cbegin(rates),
        std::find_if(std::cbegin(rates) + offset, std::cend(rates) - 1,
                     [=](const double r) { return r >= rate; })));
    return limits.at(index);
  }
}

static std::ofstream open_file(const std::string &fname) {
  if (detail::file_exists(fname)) {
    std::ostringstream os;
    os << "File " << fname << " already exists.";
    throw std::runtime_error(os.str());
  }

  std::ofstream log_file(fname.c_str(), std::ofstream::out);
  return log_file;
}

double awgn_simulation::sigma(const double eb_no) const {
  return 1.0f / sqrt((2 * decoder.rate() * pow(10, eb_no / 10.0)));
}

awgn_simulation::awgn_simulation(const class decoder &decoder_,
                                 const double step_, const uint64_t seed)
    : decoder(decoder_), step(step_), generator(seed) {}

size_t awgn_simulation::samples(const double &wer) const {
  return static_cast<size_t>(std::min(1e6, 5e3 / wer));
}

void awgn_simulation::operator()() {
  const size_t wer_width = std::numeric_limits<double>::digits10;
  const size_t ebno_width = 6;
  std::ofstream log_file(open_file(decoder.to_string() + ".log"));
  std::vector<float> b(decoder.n());

  log_file << std::setw(ebno_width + 1) << "ebno"
           << " ";
  log_file << std::setw(wer_width + 6) << "wer" << std::endl;

  const size_t tmp = static_cast<size_t>(ebno(decoder.rate()) / step);
  const double start = (tmp + (1.0 / step)) * step;
  const double max = std::max(8.0, start) + step / 2;

  double wer = 0.5;

  /* TODO round ebno() up to next step */
  for (double eb_no = start; eb_no < max; eb_no += step) {
    std::normal_distribution<float> distribution(
        1.0, static_cast<float>(sigma(eb_no)));
    auto noise = std::bind(std::ref(distribution), std::ref(generator));
    size_t word_errors = 0;
    size_t iterations = samples(wer);

    std::cout << std::this_thread::get_id() << " " << decoder.to_string()
              << ": E_b/N_0 = " << eb_no << " with " << iterations << " … ";
    std::cout.flush();
    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < iterations; i++) {
      std::generate(std::begin(b), std::end(b), noise);
      try {
        auto result = decoder.correct(b);
        if (std::any_of(std::cbegin(result), std::cend(result),
                        [](const auto &bit) { return bool(bit); })) {
          word_errors++;
        }
      }
      catch (const decoding_failure &) {
        word_errors++;
      }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(
        end_time - start_time).count();
    std::cout << seconds << " s" << std::endl;

    wer = static_cast<double>(word_errors) / iterations;

    log_file << std::setw(ebno_width + 1) << std::setprecision(ebno_width)
             << std::defaultfloat << eb_no << " ";
    log_file << std::setw(wer_width + 1) << std::setprecision(wer_width)
             << std::scientific << wer << std::endl;
  }
}

bitflip_simulation::bitflip_simulation(const class decoder &decoder_,
                                       const size_t errors_)
    : decoder(decoder_), errors(errors_) {}

void bitflip_simulation::operator()() const {
  const size_t wer_width = std::numeric_limits<double>::digits10;
  const size_t ebno_width = 6;

  std::ofstream log_file(open_file(decoder.to_string() + ".log"));

  log_file << std::setw(ebno_width + 1) << "errors"
           << " ";
  log_file << std::setw(wer_width + 6) << "wer" << std::endl;

  const size_t length = decoder.n();

  for (size_t error = 0; error <= errors; error++) {
    size_t patterns = 0;
    size_t word_errors = 0;
    std::vector<int> b;

    std::fill_n(std::back_inserter(b), length - error, 0);
    std::fill_n(std::back_inserter(b), error, 1);

    std::cout << std::this_thread::get_id() << " " << decoder.to_string()
              << ": errors = " << error << " … ";
    std::cout.flush();
    auto start = std::chrono::high_resolution_clock::now();

    do {
      std::vector<float> x;
      x.reserve(length);

      std::transform(std::cbegin(b), std::cend(b), std::back_inserter(x),
                     [](const auto &bit) { return -2 * bit + 1; });

      patterns++;
      try {
        auto result = decoder.correct(x);
        if (std::any_of(std::cbegin(result), std::cend(result),
                        [](const auto &bit) { return bool(bit); })) {
          word_errors++;
        }
      }
      catch (const decoding_failure &) {
        word_errors++;
      }
    } while (std::next_permutation(std::begin(b), std::end(b)));

    auto end = std::chrono::high_resolution_clock::now();
    auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << seconds << " s " << word_errors << " " << patterns
              << std::endl;

    log_file << std::setw(ebno_width + 1) << std::setprecision(ebno_width)
             << std::defaultfloat << error << " ";
    log_file << std::setw(wer_width + 1) << std::setprecision(wer_width)
             << std::scientific << static_cast<double>(word_errors) / patterns
             << std::endl;
  }
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"

std::atomic_bool thread_pool::running{ true };
std::mutex thread_pool::lock;
std::condition_variable thread_pool::cv;
std::list<std::function<void(void)> > thread_pool::queue;

#pragma clang diagnostic pop

thread_pool::thread_pool() {
  for (unsigned i = 0; i < std::thread::hardware_concurrency(); i++) {
    pool.emplace_back(std::bind(&thread_pool::thread_function, this));
  }
}

thread_pool::~thread_pool() {
  running = false;
  cv.notify_all();
  for (auto &&t : pool)
    t.join();
}

void thread_pool::thread_function() const {
  for (;;) {
    bool have_work = false;
    std::function<void(void)> work;
    {
      std::unique_lock<std::mutex> m(lock);

      /*
       * empty running wait
       * 0     0       0
       * 0     1       0
       * 1     0       0
       * 1     1       1
      */
      cv.wait(m, [&]() { return !(queue.empty() && running); });
      if (!queue.empty()) {
        work = queue.front();
        queue.pop_front();
        have_work = true;
      }
    }

    if (have_work) {
      try {
        work();
      }
      catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
      }
    } else if (!running) {
      return;
    }
  }
}

