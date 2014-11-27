#include <vector>
#include <stdexcept>
#include <iostream>

#include "codes/bch.h"
#include "codes/rs.h"

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (const auto &i : v)
    os << i;
  return os;
}

template <>
std::ostream &operator<<(std::ostream &os, const std::vector<uint8_t> &v) {
  for (const auto &i : v)
    os << static_cast<unsigned>(i);
  return os;
}

/* Concepts:
 * C1 operator!= C2
 * operator<< C1
 * operator<< C2
 */
template <typename C1, typename C2>
void expect_equal(const C1 &a, const C2 &b) {
  if (a != b) {
    std::cout << "Decoded word:  " << b << std::endl;
    std::cout << "Expected word: " << a << std::endl;
    throw std::runtime_error("Decoding Error.");
  }
}

static void task_6_1() {
  std::cout << std::endl << "Task 6.1" << std::endl << std::endl;

  using code_type = cyclic::primitive_bch<4, dmin<7> >;

  code_type code;

  const std::vector<unsigned char> a(
      { 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1 });

  std::vector<unsigned> b1({ 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1 });
  std::vector<unsigned> b2({ 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1 });

  std::cout << "Decoding b1" << std::endl;
  expect_equal(a, code.correct(b1));
  std::cout << "Decoding b2" << std::endl;
  expect_equal(a, code.correct(b2));
}

static void task_6_2() {
  std::cout << std::endl << "Task 6.2" << std::endl << std::endl;

  using code_type = cyclic::primitive_bch<4, dmin<5> >;

  code_type code;

  const std::vector<unsigned> a(
      { 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });
  const std::vector<unsigned> b(
      { 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0 });

  try {
    auto b_ = code.correct(b);
    std::cout << "Decoded b into: " << b << std::endl;
  }
  catch (const decoding_failure &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
    return;
  }

  throw std::runtime_error("Expected decoding failure");
}

static void task_6_3() {
  std::cout << std::endl << "Task 6.3" << std::endl << std::endl;

  using code_type = cyclic::primitive_bch<4, dmin<6> >;

  const std::vector<unsigned> a(
      { 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });

  code_type code;

  const std::vector<unsigned> b1(
      { 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1 });
  const std::vector<unsigned> b2(
      { 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1 });
  const std::vector<unsigned> b3(
      { 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });

  expect_equal(a, code.correct<unsigned>(b1));
  expect_equal(a, code.correct<unsigned>(b2));
  try {
    code.correct(b3);
  }
  catch (const std::exception &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

static void task_6_4() {
  std::cout << std::endl << "Task 6.4" << std::endl << std::endl;

  using code_type = cyclic::rs<3, errors<1> >;
  using Element = typename code_type::Element;

  std::vector<Element> b1({ Element(0),            Element(0), Element(0),
                            Element(0),            Element(0), Element(0),
                            Element::from_power(4) });

  std::vector<Element> b2({ Element::from_power(2), Element::from_power(2),
                            Element(1),             Element(0),
                            Element(0),             Element(0),
                            Element::from_power(4) });

  code_type code;

  auto b_ = code.correct(b1);
  std::cout << "Corrected vector: " << b_ << std::endl;

  try {
    b_ = code.correct(b2);
    std::cout << "Expected decoding failure." << std::endl;
  }
  catch (const decoding_failure &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

static void task_6_5() {
  std::cout << std::endl << "Task 6.5" << std::endl << std::endl;

  using RS_Code = cyclic::rs<4, errors<3> >;
  using Element = typename RS_Code::Element;

  std::vector<Element> b({ Element(1), Element(1), Element(1), Element(1),
                           Element(0), Element(0), Element(0), Element(0),
                           Element(0), Element(0), Element(0), Element(0),
                           Element(0), Element(0), Element(0) });

  RS_Code code;

  std::cout << b << std::endl;

  try {
    auto b_ = code.correct<Element>(b);
    std::cout << "Expected decoding failure." << std::endl;
  }
  catch (const decoding_failure &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

static void task_6_6() {
  std::cout << std::endl << "Task 6.6" << std::endl << std::endl;

  using RS_Code = cyclic::rs<3, errors<2> >;
  using Element = RS_Code::Element;

  std::vector<Element> b({ Element::from_power(6), Element::from_power(2),
                           Element::from_power(2), Element::from_power(5),
                           Element(0),             Element(0),
                           Element::from_power(5) });

  std::cout << b << std::endl;

  RS_Code code;

  auto b_ = code.correct<Element>(b);
  std::cout << "Corrected vector: " << b_ << std::endl;
}

static void task_6_7() {
  std::cout << std::endl << "Task 6.7" << std::endl << std::endl;

  using RS_Code = cyclic::rs<3, errors<2>, cyclic::berlekamp_massey_tag>;
  using Element = RS_Code::Element;

  const std::vector<unsigned> powers = { 6, 2, 2, 5, 4, 6, 5 };
  const std::vector<unsigned> erasures = { 5, 4, 3, 2 };

  std::vector<Element> a;
  std::transform(
      std::begin(powers), std::end(powers), std::back_inserter(a),
      [](const unsigned &power) { return Element::from_power(power); });
  std::vector<Element> b(a);

  for (const auto &erasure : erasures)
    b.at(erasure) = Element(0);

  std::cout << "a: " << a << std::endl;
  std::cout << "b: " << b << std::endl;

  RS_Code code;

  expect_equal(a, code.correct<Element>(b, erasures));
}

static void task_6_8() {
  std::cout << std::endl << "Task 6.8" << std::endl << std::endl;

  using RS_Code = cyclic::rs<3, errors<2>, cyclic::berlekamp_massey_tag>;
  using Element = RS_Code::Element;

  const std::vector<unsigned> powers = { 2, 0, 4, 0, 5, 0, 2 };
  const std::vector<unsigned> erasures = { 1, 3 };

  std::vector<Element> b;
  std::transform(
      std::begin(powers), std::end(powers), std::back_inserter(b),
      [](const unsigned &power) { return Element::from_power(power); });

  std::cout << "b: " << b << std::endl;

  RS_Code code;

  auto b_ = code.correct<Element>(b, erasures);
  std::cout << "Corrected vector: " << b_ << std::endl;
}

static void task_6_9() {
  std::cout << std::endl << "Task 6.9" << std::endl << std::endl;

  using RS_Code_pgz = cyclic::rs<3, errors<2> >;
  using RS_Code_bm = cyclic::rs<3, errors<2>, cyclic::berlekamp_massey_tag>;

  using Element = RS_Code_pgz::Element;

  std::vector<Element> b({ Element::from_power(3), Element::from_power(4),
                           Element::from_power(0), Element::from_power(3),
                           Element::from_power(4), Element::from_power(3),
                           Element::from_power(3) });

  std::cout << "b: " << b << std::endl;

  RS_Code_pgz pgz;
  RS_Code_bm bm;

  auto b_bm = bm.correct<Element>(b);
  auto b_pzg = pgz.correct<Element>(b);

  std::cout << "Corrected vector (PZG):        " << b_pzg << std::endl;
  std::cout << "Corrected vector (BMA/EUKLID): " << b_bm << std::endl;
}

static void task_6_10() {
  std::cout << std::endl << "Task 6.10" << std::endl << std::endl;

  using code_type_peterson = cyclic::primitive_bch<
      4, errors<2>, cyclic::peterson_gorenstein_zierler_tag>;
  using code_type_berlekamp =
      cyclic::primitive_bch<4, errors<2>, cyclic::berlekamp_massey_tag>;
  using code_type_euklid =
      cyclic::primitive_bch<4, errors<2>, cyclic::euklid_tag>;

  const std::vector<uint8_t> a({ 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 });

  code_type_peterson peterson;
  code_type_berlekamp berlekamp_massey;
  code_type_euklid euklid;

  std::cout << "PGZ:    " << peterson.correct(a) << std::endl;
  std::cout << "BM:     " << berlekamp_massey.correct(a) << std::endl;
  std::cout << "EUKLID: " << euklid.correct(a) << std::endl;
}

int main() {
  task_6_1();
  task_6_2();
  task_6_3();

  task_6_4();
  task_6_5();
  task_6_6();
  task_6_7();
  task_6_8();
  task_6_9();

  task_6_10();
}
