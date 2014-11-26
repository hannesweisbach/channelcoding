#pragma once

#include <vector>
#include <memory>
#include <string>
#include <utility>
#include <iterator>

#include <thread>
#include <random>
#include <functional>
#include <list>
#include <fstream>

#include <mutex>
#include <condition_variable>
#include <atomic>

/* TODO: InputSequence concept. */

class decoder {
  class decoder_concept {
  public:
    virtual ~decoder_concept() = default;
    virtual std::vector<unsigned char>
    correct(const std::vector<float> &b) const = 0;
    virtual std::string to_string() const = 0;
    virtual double rate() const = 0;
    virtual unsigned n() const = 0;
  };

  template <typename T> class decoder_model : public decoder_concept {
    T implementation;

  public:
    decoder_model(T arg) : implementation(std::move(arg)) {}
    std::vector<unsigned char> correct(const std::vector<float> &b) const
        override {
      return implementation.template correct<unsigned char>(b);
    }
    std::string to_string() const override {
      return implementation.to_string();
    }
    double rate() const override { return implementation.rate; }
    unsigned n() const override { return implementation.n; }
  };

  std::shared_ptr<const decoder_concept> _self;

public:
  template <typename T>
  decoder(T decoder)
      : _self(std::make_shared<decoder_model<T> >(std::move(decoder))) {}

  template <typename InputSequence>
  std::vector<unsigned char> correct(const InputSequence &b) const {
    return _self->correct(b);
  }
  std::string to_string() const { return _self->to_string(); }
  double rate() const { return _self->rate(); }
  unsigned n() const { return _self->n(); }
};

class awgn_simulation {
  const decoder &decoder;
  const double step;
  const uint64_t seed;
  std::mt19937_64 generator;

  double sigma(const double eb_no) const;

public:
  awgn_simulation(const class decoder &decoder, const double step = 0.5,
                  const uint64_t seed = 0);
  void operator()();
};

class bitflip_simulation {
  const decoder &decoder;
  const size_t errors = 0;

public:
  bitflip_simulation(const class decoder &decoder, const size_t errors = 0);
  void operator()() const;
};

class thread_pool {
  std::vector<std::thread> pool;

  static std::atomic_bool running;
  static std::mutex lock;
  static std::condition_variable cv;
  static std::list<std::function<void(void)> > queue;

  void thread_function() const;

public:
  thread_pool();
  ~thread_pool();

  template <typename Functor> void push(Functor &&f) {
    std::list<std::function<void(void)> > tmp;
    tmp.emplace_back(std::forward<Functor>(f));
    {
      std::lock_guard<std::mutex> m(lock);
      queue.splice(std::cend(queue), tmp);
    }
    cv.notify_one();
  }
};
