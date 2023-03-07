#pragma once

#include <iosfwd>
#include <string>
#include <vector>
#include <iostream>

struct big_integer {

public:
  big_integer();
  big_integer(big_integer const& other);
  big_integer(short a);
  big_integer(int a);
  big_integer(long a);
  big_integer(long long a);
  big_integer(unsigned short a);
  big_integer(unsigned int a);
  big_integer(unsigned long a);
  big_integer(unsigned long long a);
  explicit big_integer(std::string const& str);
  ~big_integer();

  big_integer& operator=(big_integer const& other);

  big_integer& operator+=(big_integer const& rhs);
  big_integer& operator-=(big_integer const& rhs);
  big_integer& operator*=(big_integer const& rhs);
  big_integer& operator/=(big_integer const& rhs);
  big_integer& operator%=(big_integer const& rhs);

  big_integer& operator&=(big_integer const& rhs);
  big_integer& operator|=(big_integer const& rhs);
  big_integer& operator^=(big_integer const& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend bool operator==(big_integer const& a, big_integer const& b);
  friend bool operator!=(big_integer const& a, big_integer const& b);
  friend bool operator<(big_integer const& a, big_integer const& b);
  friend bool operator>(big_integer const& a, big_integer const& b);
  friend bool operator<=(big_integer const& a, big_integer const& b);
  friend bool operator>=(big_integer const& a, big_integer const& b);

  friend std::string to_string(big_integer const& a);

private:

  static bool less_equals_abs(const big_integer& a, const big_integer& b);
  void swap(big_integer& other);
  void swap(big_integer&& other);
  void extend(size_t length_);
  void delete_leading_zeroes();
  static std::vector<uint32_t> subtract_abs(const big_integer& bigger,
                                            const big_integer& smaller);
  void plus_minus(big_integer const& rhs, bool flag);
  void abs_increase(big_integer const& rhs);
  static big_integer product(big_integer const& y, uint32_t k);
  void mul_uint32_t(uint32_t k);
  void abs_increase_uint32_t(uint32_t k);
  void abs_decrease_uint32_t(uint32_t k);
  void quotient(uint32_t k);
  uint32_t remainder(uint32_t k);
  uint32_t trial(big_integer const& d);
  bool smaller(big_integer const& dq, uint32_t m);
  void difference(big_integer const& dq, uint32_t m);
  std::pair<big_integer, big_integer> long_divide(big_integer const& y);

  void to_form_with_addition();

  void binary_op(big_integer const& rhs,
                 uint32_t (*func)(uint32_t& a, uint32_t& b), bool tmp_sign);

  std::vector<uint32_t> data_;
  bool sign = false;

  // for of store number
  // number : sign = number < 0,
  // number_abs = sum i = 0 .. data.size data[i] * (base ^ i)
  // number = 0 : sign = false && data_.empty()
};

big_integer operator+(big_integer a, big_integer const& b);
big_integer operator-(big_integer a, big_integer const& b);
big_integer operator*(big_integer a, big_integer const& b);
big_integer operator/(big_integer a, big_integer const& b);
big_integer operator%(big_integer a, big_integer const& b);

big_integer operator&(big_integer a, big_integer const& b);
big_integer operator|(big_integer a, big_integer const& b);
big_integer operator^(big_integer a, big_integer const& b);

big_integer operator<<(big_integer a, int b);
big_integer operator>>(big_integer a, int b);

bool operator==(big_integer const& a, big_integer const& b);
bool operator!=(big_integer const& a, big_integer const& b);
bool operator<(big_integer const& a, big_integer const& b);
bool operator>(big_integer const& a, big_integer const& b);
bool operator<=(big_integer const& a, big_integer const& b);
bool operator>=(big_integer const& a, big_integer const& b);

std::string to_string(big_integer const& a);
std::ostream& operator<<(std::ostream& s, big_integer const& a);