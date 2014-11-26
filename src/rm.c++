#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iterator>

void print_matrix(const std::vector<std::vector<int> > &matrix) {
  for (const auto &row : matrix) {
    for (const auto &e : row)
      std::cout << e << " ";
    std::cout << std::endl;
  }
}

std::vector<int> wedge_product(const std::vector<int> &lhs,
                               const std::vector<int> &rhs) {
  if (rhs.size() != lhs.size()) {
    std::ostringstream os;
    os << "Vectors have different sizes of " << lhs.size() << " and "
       << rhs.size();
    throw std::runtime_error(os.str());
  }

  std::vector<int> result;

  auto it1 = std::begin(lhs);
  auto it2 = std::begin(rhs);

  std::generate_n(std::back_inserter(result), lhs.size(), [&]() {
    auto v1 = *it1;
    auto v2 = *it2;
    ++it1;
    ++it2;
    return v1 * v2;
  });

  return result;
}

int main() {
  const unsigned order = 1;
  const unsigned m = 13;

  const unsigned length = 2 << (m - 1);
  std::vector<std::vector<int> > matrix;

  matrix.push_back(std::vector<int>(length, 1));

  for (unsigned i = 0; i < m; i++) {
    std::vector<int> v;
    const unsigned repeat = 1 << i;
    const unsigned blocks = length / repeat;
    unsigned value = 1;
    for (unsigned block = 0; block < blocks; block++) {
      for (unsigned r = 0; r < repeat; r++) {
        v.push_back(value);
      }
      value ^= 1;
    }
    matrix.push_back(v);
  }

  print_matrix(matrix);
  std::vector<std::vector<int> > tmp;

  for (auto first_start = std::cbegin(matrix) + 1;
       first_start != std::cend(matrix) - 1; ++first_start) {
    for (auto next = first_start + 1; next != std::cend(matrix); ++next) {
      tmp.push_back(wedge_product(*first_start, *next));
    }
  }

  matrix.insert(std::end(matrix), std::begin(tmp), std::end(tmp));
  print_matrix(matrix);
}
