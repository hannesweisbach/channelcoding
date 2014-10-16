#include <stdexcept>
#include <iostream>

#include "gf.h"
#include "bch.h"
#include "rs.h"
#include "util.h"

void expect_equal(const std::vector<int> &a, const std::vector<int> &b) {
  if (!std::equal(std::begin(b), std::end(b), std::begin(a), std::end(a))) {
    std::cout << "Decoded word:   ";
    for (const auto &i : b)
      std::cout << i << " ";
    std::cout << std::endl;
    std::cout << "Expected word: ";
    for (const auto &i : a)
      std::cout << i << " ";
    throw std::runtime_error("Decoding Error.");
  }
}

void expect_equal(const gf_polynomial &a, const gf_polynomial &b) {
  if (a != b) {
    std::cout << "Decoded word:   " << b << std::endl;
    std::cout << "Expected word: " << a << std::endl;
    throw std::runtime_error("Decoding Error.");
  }
}

void task_6_1() {
  std::cout << std::endl << "Task 6.1" << std::endl << std::endl;
  
  const std::vector<int> a({ 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1 });
  bch code(4, 0x13, 7);
  std::vector<int> b1({ 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1 });
  std::vector<int> b2({ 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1 });

  std::cout << "Decoding b1" << std::endl;
  expect_equal(a, code.correct_peterson(b1));
  std::cout << "Decoding b2" << std::endl;
  expect_equal(a, code.correct_peterson(b2));
}

void task_6_2() {
  std::cout << std::endl << "Task 6.2" << std::endl << std::endl;
  
  bch code(4, 0x13, 5);
  const std::vector<int> a({ 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });
  std::vector<int> b({ 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0 });

  try {
    auto b_ = code.correct_peterson(b);
    std::cout << "Decoded b into: ";
    for (const auto &i : b_)
      std::cout << i << " ";
  }
  catch (const std::exception &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
    return;
  }

  throw std::runtime_error("Expected decoding failure");
}

void task_6_3() {
  std::cout << std::endl << "Task 6.3" << std::endl << std::endl;

  const std::vector<int> a({ 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });
  bch code(4, 0x13, 6);
  std::vector<int> b1({ 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1 });
  std::vector<int> b2({ 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1 });
  std::vector<int> b3({ 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1 });

  expect_equal(a, code.correct_peterson(b1));
  expect_equal(a, code.correct_peterson(b2));
  try {
    code.correct_peterson(b3);
  }
  catch (const std::exception &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

void task_6_4() {
  std::cout << std::endl << "Task 6.4" << std::endl << std::endl;
  
  gf field(3, 0xb);
  
  gf_polynomial b1(field, std::vector<gf_element>(7, field.zero()));
  b1.at(6) = field.power_to_polynomial(4);
  
  rs code(3, 0xb, 3);
  
  auto b_ = code.correct_pzg(b1);
  std::cout << "Corrected vector: " << b_ << std::endl;

  gf_polynomial b2(field, std::vector<gf_element>(7, field.zero()));
  b2.at(1) = field.power_to_polynomial(2);
  b2.at(6) = field.power_to_polynomial(4);

  try {
    b_ = code.correct_pzg(b2);
    std::cout << "Expected decoding failure." << std::endl;
  }
  catch (const std::exception &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

void task_6_5() {
  std::cout << std::endl << "Task 6.5" << std::endl << std::endl;
  
  gf field(4, 0x13);
  
  gf_polynomial b(field, std::vector<gf_element>(15, field.zero()));
  for(int i = 0; i < 4; i++)
    b.at(i) = field.power_to_polynomial(0);
  
  std::cout << b << std::endl;
  rs code(4, 0x13, 7);

  try {
    auto b_ = code.correct_pzg(b);
    std::cout << "Expected decoding failure." << std::endl;
  }
  catch (const std::exception &e) {
    std::cout << "Decoding failure:" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

void task_6_6() {
  std::cout << std::endl << "Task 6.6" << std::endl << std::endl;

  gf field(3, 0xb);

  gf_polynomial b(field, std::vector<gf_element>(7, field.zero()));
  b.at(6) = field.power_to_polynomial(5);
  b.at(3) = field.power_to_polynomial(5);
  b.at(2) = field.power_to_polynomial(2);
  b.at(1) = field.power_to_polynomial(2);
  b.at(0) = field.power_to_polynomial(6);

  std::cout << b << std::endl;
  rs code(3, 0xb, 5);

  auto b_ = code.correct_pzg(b);
  std::cout << "Corrected vector: " << b_ << std::endl;
}

void task_6_7() {
  std::vector<int> powers = { 6, 2, 2, 5, 4, 6, 5 };
  std::vector<int> erasures = { 5, 4, 3, 2 };
  std::cout << std::endl << "Task 6.7" << std::endl << std::endl;

  gf field(3, 0xb);
  
  gf_polynomial a(field, powers, true);
  gf_polynomial b(a);

  for (const auto &erasure : erasures)
    b.at(erasure) = field.one();

  std::cout << "a: " << a << std::endl;
  std::cout << "b: " << b << std::endl;

  rs code(3, 0xb, 5);

  expect_equal(a,
               code.correct_erasures(b, gf_polynomial(field, erasures, true)));
}

void task_6_8() {
  std::vector<int> powers = { 2, 0, 4, 0, 5, 0, 2 };
  std::vector<int> erasures = { 1, 3 };
  std::cout << std::endl << "Task 6.8" << std::endl << std::endl;

  gf field(3, 0xb);

  gf_polynomial b(field, powers, true);
  b.at(5) = field.zero();

  for (const auto &erasure : erasures)
    b.at(erasure) = field.one();

  std::cout << "b: " << b << std::endl;

  rs code(3, 0xb, 5);

  auto b_ = code.correct_bm(b, gf_polynomial(field, erasures, true));
  std::cout << "Corrected vector: " << b_ << std::endl;
}

void task_6_9() {
  std::vector<int> powers = { 3, 4, 0, 3, 4, 3, 3 };
  std::cout << std::endl << "Task 6.9" << std::endl << std::endl;

  gf field(3, 0xb);

  gf_polynomial b(field, powers, true);

  std::cout << "b: " << b << std::endl;

  rs code(3, 0xb, 5);

  auto b_bm = code.correct_bm(b);
  auto b_pzg = code.correct_pzg(b);

  std::cout << "Corrected vector (PZG):        " << b_pzg << std::endl;
  std::cout << "Corrected vector (BMA/EUKLID): " << b_bm << std::endl;
}

void task_6_10() {
  std::cout << std::endl << "Task 6.10" << std::endl << std::endl;

  const std::vector<int> coefficients(
      { 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 });
  bch code(4, 0x13, 5);
  
  gf field(4, 0x13);
  gf_polynomial b(field, coefficients);

  //std::cout << "Corrected: " <<;
  code.correct_peterson(coefficients);
  code.correct_bm(b);
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

#if 0
  std::vector<int> powers = { 6, 2, 3, 5, 4, 6, 5 };
  std::vector<gf_element> erasures;

  gf field(3, 0xb);
  gf_polynomial a(field, powers, true);
  gf_polynomial b(a);

  std::cout << "a: " << a << std::endl;
  std::cout << "b: " << b << std::endl;

  rs code(3, 0xb, 5);
  
  expect_equal(a, code.correct_pzg(b));

  //expect_equal(a, code.correct_erasures(b, erasures));
#endif
}
