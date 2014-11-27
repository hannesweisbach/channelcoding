#pragma once

#include <cstdlib>
#include <iomanip>

template <typename charT, typename traits = std::char_traits<charT> >
class center_helper {
  const std::basic_string<charT, traits> &str_;

public:
  center_helper(const std::basic_string<charT, traits> &str) : str_(str) {}
  template <typename a, typename b>
  friend std::basic_ostream<a, b> &operator<<(std::basic_ostream<a, b> &s,
                                              const center_helper<a, b> &c);
};

template <typename charT, typename traits = std::char_traits<charT> >
center_helper<charT, traits>
centered(const std::basic_string<charT, traits> &str) {
  return center_helper<charT, traits>(str);
}

// redeclare for std::string directly so we can support anything that implicitly
// converts to std::string
inline center_helper<std::string::value_type, std::string::traits_type>
centered(const std::string &str) {
  return center_helper<char>(str);
}

template <typename charT, typename traits = std::char_traits<charT> >
std::basic_ostream<charT, traits> &operator<<(
    std::basic_ostream<charT, traits> &s,
    const center_helper<charT, traits> &c) {
  std::streamsize w = s.width();
  ssize_t strlen = mbstowcs(nullptr, c.str_.c_str(), 0);
  if (w > strlen) {
    std::streamsize left = (w + strlen) / 2;
    s << std::setw(left - strlen) << "" << c.str_ << std::setw(w - left) << "";
  } else {
    s << c.str_;
  }
  return s;
}


