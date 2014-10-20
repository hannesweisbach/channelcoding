#pragma once

#ifdef __linux__
namespace std {
template <typename T> auto rbegin(T &c) -> decltype(c.rbegin()) {
  return c.rbegin();
}
template <typename T> auto rend(T &c) -> decltype(c.rend()) { return c.rend(); }

template <typename T> auto cbegin(const T &c) -> decltype(c.cbegin()) {
  return c.cbegin();
}
template <typename T> auto crbegin(const T &c) -> decltype(c.crbegin()) {
  return c.crbegin();
}
template <typename T> auto cend(const T &c) -> decltype(c.cend()) {
  return c.cend();
}
template <typename T> auto crend(const T &c) -> decltype(c.crend()) {
  return c.crend();
}
}
#endif

#include <string>
#include <sys/stat.h>

class decoding_failure : std::runtime_error {
public:
  decoding_failure(const std::string &what) : runtime_error(what) {}
};

inline bool file_exists(const std::string &fname) {
  struct stat buf;
  return (stat(fname.c_str(), &buf) != -1);
}
