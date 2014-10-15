#pragma once

namespace std {
template <typename T> auto cbegin(const T &c) { return c.cbegin(); }
template <typename T> auto crbegin(const T &c) { return c.crbegin(); }
template <typename T> auto cend(const T &c) { return c.cend(); }
template <typename T> auto crend(const T &c) { return c.crend(); }
}

