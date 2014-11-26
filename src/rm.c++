#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iterator>

void print_matrix(const std::vector<std::vector<int> > &matrix);
std::vector<int> wedge_product(const std::vector<int> &lhs,
                               const std::vector<int> &rhs);

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
  result.reserve(rhs.size());

  std::transform(std::cbegin(lhs), std::cend(lhs), std::begin(rhs),
                 std::back_inserter(result), std::multiplies<int>{});

  return result;
}

int main() {
  // const unsigned order = 1;
  const unsigned m = 13;

  const unsigned length = 2 << (m - 1);
  std::vector<std::vector<int> > matrix;

  /* 1st order */
  matrix.push_back(std::vector<int>(length, 1));

  for (unsigned i = 0; i < m; i++) {
    std::vector<int> v;
    const unsigned repeat = 1 << i;
    const unsigned blocks = length / repeat;
    int value = 1;
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

  /* 2nd order */
  for (auto first_start = std::cbegin(matrix) + 1;
       first_start != std::cend(matrix) - 1; ++first_start) {
    for (auto next = first_start + 1; next != std::cend(matrix); ++next) {
      tmp.push_back(wedge_product(*first_start, *next));
    }
  }

  matrix.insert(std::end(matrix), std::begin(tmp), std::end(tmp));
  print_matrix(matrix);
}
