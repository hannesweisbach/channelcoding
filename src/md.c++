#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <tuple>
#include <vector>
#include <numeric>
#include <iterator>
#include <string>

void print_sequence(const std::vector<long> &sequence);
long euklidean_distance(const std::vector<long> &x, const std::vector<long> &y);

void print_sequence(const std::vector<long> &sequence) {
  std::cout << "(";
  for (auto it = std::cbegin(sequence); it != std::cend(sequence) - 1; ++it) {
    std::cout << *it << ", ";
  }
  std::cout << sequence.back() << ")";
}

long euklidean_distance(const std::vector<long> &x,
                        const std::vector<long> &y) {
  if (x.size() != y.size())
    throw std::runtime_error("Sizes do not match");
  return std::inner_product(
      std::begin(x), std::end(x), std::begin(y), 0L,
      [](const auto &lhs, const auto &rhs) { return rhs + lhs; },
      [](const auto &lhs,
         const auto &rhs) { return (lhs - rhs) * (lhs - rhs); });
}

int main(int argc, const char *const argv[]) {
  const std::vector<std::vector<long> > values(
      { { 3, 3 }, { 3, -4 }, { -4, 3 }, { -4, -4 } });
  for (int arg = 1; arg + 1 < argc; arg += 2) {
    std::vector<long> bit;
    auto first = std::stol(argv[arg]);
    auto second = std::stol(argv[arg + 1]);
    bit.push_back(first);
    bit.push_back(second);
    std::cout << "Received ";
    print_sequence(bit);
    std::cout << std::endl;
    for (const auto &value : values) {
      std::cout << "Euklidean distance to ";
      print_sequence(value);
      std::cout << " = " << euklidean_distance(value, bit) << std::endl;
    }
  }
}

