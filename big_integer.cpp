#include "big_integer.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <cstddef>
#include <stdexcept>

const uint64_t BASE = static_cast<uint64_t>(1) << 32;

static uint32_t and_(uint32_t& a, uint32_t& b) {
  return a & b;
}

static uint32_t or_(uint32_t& a, uint32_t& b) {
  return a | b;
}

static uint32_t xor_(uint32_t& a, uint32_t& b) {
  return a ^ b;
}

big_integer::big_integer() = default;

big_integer::big_integer(big_integer const& other) = default;

big_integer::big_integer(short a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(int a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(long a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(long long a) : sign(a < 0) {
  bool flag = false;
  if (a == std::numeric_limits<long long>::min()) {
    a++;
    flag = true;
  }
  if (a != 0) {
    a = a < 0 ? -a : a;
    data_.push_back(static_cast<uint32_t>(a & UINT32_MAX));
    if (a > UINT32_MAX) {
      data_.push_back(static_cast<uint32_t>(a >> 32));
    }
    if (flag) {
      (*this)--;
    }
  }
}

big_integer::big_integer(unsigned short a)
    : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(unsigned int a)
    : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(unsigned long a)
    : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(unsigned long long a) {
  if (a > 0) {
    data_.push_back(static_cast<uint32_t>(a & UINT32_MAX));
    if (a > UINT32_MAX) {
      data_.push_back(static_cast<uint32_t>(a >> 32));
    }
  }
}

big_integer::big_integer(std::string const& str) : big_integer() {
  if (str.empty() || str == "-" || str == "+") {
    throw std::invalid_argument(str);
  }
  for (size_t i = str[0] == '-' || str[0] == '+' ? 1 : 0; i < str.length();) {
    size_t j = 0;
    uint32_t to_add = 0;
    uint32_t to_mul = 1;
    // read block 1e9 < uint32_max
    while (j < 9 && i + j < str.length()) {
      if (str[i + j] < '0' || str[i + j] > '9') {
        throw std::invalid_argument(str);
      }
      to_add *= 10;
      to_add += static_cast<uint32_t>(str[i + j] - '0');
      j++;
      to_mul *= 10;
    }
    i += j;
    mul_uint32_t(to_mul);
    abs_increase_uint32_t(to_add);
  }
  sign = str[0] == '-';
  delete_leading_zeroes();
}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(big_integer const& other) {
  if (this != &other) {
    big_integer(other).swap(*this);
  }
  return *this;
}

std::vector<uint32_t> big_integer::subtract_abs(const big_integer& bigger,
                                                const big_integer& smaller) {
  std::vector<uint32_t> res(bigger.data_);
  bool carry = false;
  for (size_t i = 0; i != res.size(); ++i) {
    if (carry) {
      carry = res[i] == 0;
      res[i]--;
    }
    if (i < smaller.data_.size()) {
      if (res[i] < smaller.data_[i]) {
        carry = true;
      }
      res[i] -= smaller.data_[i];
    }
  }
  return res;
}

void big_integer::plus_minus(big_integer const& rhs, bool flag) {
  if (flag) {
    abs_increase(rhs);
  } else {
    // subtract
    if (less_equals_abs(*this, rhs)) {
      sign = !sign;
      data_ = subtract_abs(rhs, *this);
    } else {
      // sign = sign;
      data_ = subtract_abs(*this, rhs);
    }
  }
  delete_leading_zeroes();
}

big_integer& big_integer::operator+=(big_integer const& rhs) {
  plus_minus(rhs, sign == rhs.sign);
  return *this;
}

big_integer& big_integer::operator-=(big_integer const& rhs) {
  plus_minus(rhs, sign != rhs.sign);
  return *this;
}

big_integer& big_integer::operator*=(big_integer const& rhs) {
  if (data_.empty() || rhs.data_.empty()) {
    return *this = 0;
  }
  big_integer res;
  res.sign = (sign != rhs.sign);
  res.extend(data_.size() + rhs.data_.size());
  uint64_t carry;
  // like in wiki
  for (size_t i = 0; i != data_.size(); ++i) {
    carry = 0;
    for (size_t j = 0; j != rhs.data_.size(); ++j) {
      uint64_t a = static_cast<uint64_t>(data_[i]) * rhs.data_[j]
                 + carry + res.data_[i + j];
      carry = a >> 32;
      res.data_[i + j] = a & UINT32_MAX;
    }
    res.data_[i + rhs.data_.size()] += carry;
  }
  (*this) = res;
  delete_leading_zeroes();
  return *this;
}

big_integer big_integer::product(const big_integer& y, uint32_t k) {
  big_integer x;
  x.extend(y.data_.size());
  uint32_t carry = 0;
  for (size_t i = 0; i < y.data_.size(); ++i) {
    uint64_t tmp = static_cast<uint64_t>(y.data_[i]) * k + carry;
    x.data_[i] = tmp % BASE;
    carry = tmp / BASE;
  }
  if (carry != 0) {
    x.data_.push_back(carry);
  }
  return x;
}

void big_integer::mul_uint32_t(uint32_t k) {
  uint32_t carry = 0;
  for (size_t i = 0; i < data_.size(); ++i) {
    uint64_t tmp = static_cast<uint64_t>(data_[i]) * k + carry;
    data_[i] = tmp % BASE;
    carry = tmp / BASE;
  }
  if (carry != 0) {
    data_.push_back(carry);
  }
}

void big_integer::quotient(uint32_t k) {
  uint32_t carry = 0;
  for (size_t i = data_.size(); i > 0; --i) {
    uint64_t tmp = carry * BASE + data_[i - 1];
    data_[i - 1] = tmp / k;
    carry = tmp % k;
  }
  delete_leading_zeroes();
}

uint32_t big_integer::remainder(uint32_t k) {
  uint32_t carry = 0;
  for (size_t i = data_.size(); i > 0; --i) {
    carry = (carry * BASE + data_[i - 1]) % k;
  }
  return carry;
}

uint32_t big_integer::trial(const big_integer& d) {
  uint64_t r3 = (data_.back() * BASE + data_[data_.size() - 2]);
  uint64_t d2 = d.data_.back();
  return std::min(static_cast<uint32_t>(r3 / d2), UINT32_MAX);
}

bool big_integer::smaller(const big_integer& dq, uint32_t m) {
  uint32_t i = 1;
  while (i <= data_.size()) {
    uint32_t a = (m - i < dq.data_.size() ? dq.data_[m - i] : 0);
    if (data_[data_.size() - i] != a) {
      return data_[data_.size() - i] < a;
    }
    i++;
  }
  return false;
}

void big_integer::difference(const big_integer& dq, uint32_t m) {
  uint32_t borrow = 0;
  uint32_t k = data_.size() - m;
  for (size_t i = 0; i != m; ++i) {
    uint64_t diff = BASE + data_[i + k] - dq.data_[i] - borrow;
    data_[i + k] = diff % BASE;
    borrow = 1 - diff / BASE;
  }
}

std::pair<big_integer, big_integer> big_integer::long_divide(big_integer const& y) {
  bool q_sign = sign ^ y.sign;
  bool sgn = sign;
  if (less_equals_abs(*this, y)) {
    return {0, *this};
  }
  big_integer q;
  uint32_t f = BASE / (y.data_[y.data_.size() - 1] + 1);
  *this *= f;
  data_.push_back(0);
  q.extend(data_.size() - y.data_.size());
  big_integer d = product(y, f);
  for (size_t k = q.data_.size(); k > 0; --k) {
    uint32_t qt = trial(d);
    big_integer dq = product(d, qt);
    while (smaller(dq, y.data_.size() + 1)) {
      qt--;
      dq -= d;
    }
    q.data_[k - 1] = qt;
    dq.data_.push_back(0);
    difference(dq, y.data_.size() + 1);
    dq.data_.pop_back();
    data_.pop_back();
  }
  quotient(f);
  q.sign = q_sign;
  sign = sgn;
  q.delete_leading_zeroes();
  delete_leading_zeroes();
  return {q, *this};
}

big_integer& big_integer::operator/=(big_integer const& rhs) {
  if (rhs.data_.empty()) {
    throw std::invalid_argument("division by zero");
  }
  swap(long_divide(rhs).first);
  return *this;
}

big_integer& big_integer::operator%=(big_integer const& rhs) {
  if (rhs.data_.empty()) {
    throw std::invalid_argument("division by zero");
  }
  swap(long_divide(rhs).second);
  return *this;
}

void big_integer::to_form_with_addition() {
  // pre: *this < 0
  // if number >= 0 addition form == our form
  for (uint32_t & i : data_) {
    i = UINT32_MAX - i;
  }
  sign = false;
  *this += 1;
}

void big_integer::binary_op(const big_integer& rhs,
                            uint32_t (*func)(uint32_t&, uint32_t&),
                            bool tmp_sign) {
  // extract common code
  big_integer a = rhs;
  extend(a.data_.size());
  a.extend(data_.size());
  if (sign) {
    to_form_with_addition();
  }
  if (a.sign) {
    a.to_form_with_addition();
  }
  for (size_t i = 0; i != data_.size(); ++i) {
    data_[i] = func(data_[i], a.data_[i]);
  }
  if (tmp_sign) {
    to_form_with_addition();
  }
  sign = tmp_sign;
  delete_leading_zeroes();
}

big_integer& big_integer::operator&=(big_integer const& rhs) {
  binary_op(rhs, &and_, sign & rhs.sign);
  return *this;
}

big_integer& big_integer::operator|=(big_integer const& rhs) {
  binary_op(rhs, &or_, sign | rhs.sign);
  return *this;
}

big_integer& big_integer::operator^=(big_integer const& rhs) {
  binary_op(rhs, &xor_, sign ^ rhs.sign);
  return *this;
}

big_integer& big_integer::operator<<=(int rhs) {
  if (rhs == 0) {
    return *this;
  } else if (rhs > 0) {
    mul_uint32_t(static_cast<uint32_t>(1) << (rhs % 32));
    size_t to_push = rhs / 32;
    data_.insert(data_.cbegin(), to_push, 0);
  } else {
    rhs = -rhs;
    quotient(static_cast<uint32_t>(1) << (rhs % 32));
    size_t to_delete = rhs / 32;
    if (to_delete > data_.size()) {
      to_delete = data_.size();
    }
    data_.erase(data_.begin(), data_.begin() + to_delete);
    delete_leading_zeroes();
    if (sign) {
      (*this)--;
    }
  }
  delete_leading_zeroes();
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  return (*this) <<= (-rhs);
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  big_integer neg(*this);
  if (!neg.data_.empty()) {
    neg.sign = !sign;
  }
  return neg;
}

big_integer big_integer::operator~() const {
  return -(*this) - 1;
}

big_integer& big_integer::operator++() {
  sign ? abs_decrease_uint32_t(1) : abs_increase_uint32_t(1);
  return *this;
}

big_integer big_integer::operator++(int) {
  big_integer r = *this;
  ++* this;
  return r;
}

big_integer& big_integer::operator--() {
  sign ? abs_increase_uint32_t(1) : abs_decrease_uint32_t(1);
  return *this;
}

big_integer big_integer::operator--(int) {
  big_integer r = *this;
  --* this;
  return r;
}

big_integer operator+(big_integer a, big_integer const& b) {
  return a += b;
}

big_integer operator-(big_integer a, big_integer const& b) {
  return a -= b;
}

big_integer operator*(big_integer a, big_integer const& b) {
  return a *= b;
}

big_integer operator/(big_integer a, big_integer const& b) {
  return a /= b;
}

big_integer operator%(big_integer a, big_integer const& b) {
  return a %= b;
}

big_integer operator&(big_integer a, big_integer const& b) {
  return a &= b;
}

big_integer operator|(big_integer a, big_integer const& b) {
  return a |= b;
}

big_integer operator^(big_integer a, big_integer const& b) {
  return a ^= b;
}

big_integer operator<<(big_integer a, int b) {
  return a <<= b;
}

big_integer operator>>(big_integer a, int b) {
  return a >>= b;
}

bool operator==(big_integer const& a, big_integer const& b) {
  return (a.sign == b.sign && a.data_ == b.data_);
}

bool operator!=(big_integer const& a, big_integer const& b) {
  return !(a == b);
}

bool operator<(big_integer const& a, big_integer const& b) {
  if (a.sign && !b.sign) {
    return true;
  } else if (!a.sign && b.sign) {
    return false;
  } else if (!a.sign) {
    // both non-negative
    return big_integer::less_equals_abs(a, b) && b != a;
  } else {
    // both negative
    return big_integer::less_equals_abs(b, a) && b != a;
  }
}

bool operator>(big_integer const& a, big_integer const& b) {
  return b < a;
}

bool operator<=(big_integer const& a, big_integer const& b) {
  if (a.sign && !b.sign) {
    return true;
  } else if (!a.sign && b.sign) {
    return false;
  } else if (!a.sign) {
    // both non-negative
    return big_integer::less_equals_abs(a, b);
  } else {
    // both negative
    return big_integer::less_equals_abs(b, a);
  }
}

bool operator>=(big_integer const& a, big_integer const& b) {
    return b <= a;
}

std::string to_string(big_integer const& a) {
  if (a.data_.empty()) {
    return "0";
  }
  big_integer tmp_big = a;
  std::string res;
  while (!tmp_big.data_.empty()) {
    uint32_t remainder = tmp_big.remainder(1000000000);
    std::string tmp = std::to_string(remainder);
    std::reverse(tmp.begin(), tmp.end());
    size_t to_push = 9 - tmp.length();
    tmp += std::string(to_push, '0');
    res += tmp;
    tmp_big.quotient(1000000000);
  }
  std::reverse(res.begin(), res.end());
  size_t to_delete = 0;
  while (res[to_delete] == '0' && to_delete < res.length()) {
    to_delete++;
  }
  res.erase(res.begin(), res.begin() + to_delete);
  if (a.sign) {
    res.insert(res.begin(), '-');
  }
  return res;
}

std::ostream& operator<<(std::ostream& s, big_integer const& a) {
  return s << to_string(a);
}

void big_integer::abs_increase_uint32_t(uint32_t k) {
  if (data_.empty()) {
    data_.push_back(k);
    return;
  }
  bool carry = data_[0] + k < data_[0];
  data_[0] += k;
  for (size_t i = 1; i < data_.size(); ++i) {
    if (carry) {
      if (data_[i] == std::numeric_limits<uint32_t>::max()) {
        data_[i] = 0;
      } else {
        data_[i]++;
        carry = false;
      }
    } else {
      break;
    }
  }
  if (carry) {
    data_.push_back(1);
  }
}

void big_integer::abs_decrease_uint32_t(uint32_t k) {
  if (data_[0] >= k) {
    data_[0] -= k;
  } else {
    data_[0] = k - data_[0];
    bool carry = true;
    size_t j = 1;
    while (j < data_.size() && carry) {
      if (data_[j] == 0) {
        data_[j] = std::numeric_limits<uint32_t>::max();
      } else {
        data_[j]--;
        carry = false;
      }
      j++;
    }
  }
  delete_leading_zeroes();
}

void big_integer::abs_increase(const big_integer& rhs) {
  uint64_t carry = 0;
  extend(rhs.data_.size());
  for (size_t i = 0; i < data_.size(); i++) {
    uint32_t prev_data_i = data_[i];
    data_[i] = data_[i] + (i < rhs.data_.size() ? rhs.data_[i] : 0) + carry;
    carry = data_[i] < prev_data_i ? 1 : 0;
  }
  if (carry > 0) {
    data_.push_back(carry);
  }
}

bool big_integer::less_equals_abs(const big_integer& a, const big_integer& b) {
  if (a.data_.size() < b.data_.size()) {
    return true;
  }
  if (a.data_.size() > b.data_.size()) {
    return false;
  }
  for (size_t i = a.data_.size(); i > 0; --i) {
    if (a.data_[i - 1] < b.data_[i - 1]) {
      return true;
    }
    else if (a.data_[i - 1] > b.data_[i - 1]) {
      return false;
    }
  }
  return true;
}

void big_integer::swap(big_integer& other) {
  std::swap(data_, other.data_);
  std::swap(sign, other.sign);
}

void big_integer::swap(big_integer&& other) {
  std::swap(data_, other.data_);
  std::swap(sign, other.sign);
}

void big_integer::extend(size_t length_) {
  if (data_.size() < length_) {
    size_t ext = length_ - data_.size();
    for (size_t i = 0; i != ext; ++i) {
      data_.push_back(0);
    }
  }
}

void big_integer::delete_leading_zeroes() {
  size_t i = data_.size();
  while (i > 0 && data_[i - 1] == 0) {
    data_.pop_back();
    i--;
  }
  if (data_.empty()) {
    sign = false;
  }
}
