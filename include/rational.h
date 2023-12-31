#ifndef RATIONAL_H
#define RATIONAL_H

#include "config.h"
#include <iostream>
#include <limits.h>
#include <stdlib.h>

#ifndef LONG_BITS
#define LONG_BITS ilog(LONG_MAX)
#endif

#ifndef RATIONAL_PRINT_MAX_DECIMALS
#define RATIONAL_PRINT_MAX_DECIMALS 9
#endif

BEGIN_HSPS_NAMESPACE

long euclid(long n, long k, long& a, long& b);
long gcd(long n, long k);
long lcm(long n, long k);
unsigned long ilog(unsigned long n);
long imag(long n);

class rational {
  long nm;
  long dv;

 public:
  rational();
  rational(int n);
  rational(long n);
  rational(long n, long d);
  rational(double v);
  rational(const rational& r);

  struct XR {
    long x_nm;
    long x_dv;

    XR& operator=(const rational r);
  };
  rational(const XR& x);

  long numerator() const;
  long divisor() const;
  long sign() const;
  bool zero() const;
  bool finite() const;
  bool infinite() const;
  bool integral() const;

  static rational reduce(rational r);
  static void split(const rational r, long& s, long& k, long& m);
  static rational invert(const rational r);
  static rational infinity(const rational r);
  static rational infinity(const long s);
  static rational floor(const rational r);
  static rational floor_to(const rational r, long div);
  static rational ceil(const rational r);
  static rational ceil_to(const rational r, long div);
  static rational frac(const rational r);
  static rational round(const rational r, long div_max);
  static rational min(const rational r0, const rational r1);
  static rational max(const rational r0, const rational r1);
  static rational rgcd(const rational r0, const rational r1);
  static rational dtor(double v);
  static rational ator(const char* s);

  rational reduce() const;
  void split(long& s, long& k, long& m) const;
  rational invert() const;
  rational floor() const;
  rational floor_to(long d) const;
  rational ceil() const;
  rational frac() const;
  rational round(long div_max) const;
  rational round() const;

  rational operator=(const rational r);
  rational operator=(long n);

  rational operator+=(const rational r);
  rational operator-=(const rational r);
  rational operator*=(const rational r);
  rational operator/=(const rational r);
  rational operator+=(long n);
  rational operator-=(long n);
  rational operator*=(long n);
  rational operator/=(long n);

  double decimal() const;

  struct printdecimal {
    const rational& r;
    printdecimal(const rational& _r) : r(_r) { };
  };

};

bool operator==(const rational r0, const rational r1);
bool operator==(const rational r0, long n1);
bool operator==(long n0, const rational r1);

bool operator!=(const rational r0, const rational r1);
bool operator!=(const rational r0, long n1);
bool operator!=(long n0, const rational r1);

bool operator<(const rational r0, const rational r1);
bool operator<=(const rational r0, const rational r1);
bool operator>(const rational r0, const rational r1);
bool operator>=(const rational r0, const rational r1);
rational operator+(const rational r0, const rational r1);
rational operator-(const rational r0, const rational r1);
rational operator*(const rational r0, const rational r1);
rational operator/(const rational r0, const rational r1);
rational operator+(const rational r0, long n1);
rational operator-(const rational r0, long n1);
rational operator*(const rational r0, long n1);
rational operator/(const rational r0, long n1);
rational operator+(long n0, const rational r1);
rational operator-(long n0, const rational r1);
rational operator*(long n0, const rational r1);
rational operator/(long n0, const rational r1);

rational safeadd(const rational r0, const rational r1);
rational safemul(const rational r0, const rational r1);

// default print op prints the rational in rational form
::std::ostream& operator<<(::std::ostream& s, const rational r);

// prints the rational in decimal form
::std::ostream& operator<<(::std::ostream& s, const rational::printdecimal& r);

// inlines

#include "config.h"

inline rational::XR& rational::XR::operator=(const rational r)
{
  x_nm = r.numerator();
  x_dv = r.divisor();
  return *this;
}

inline rational::rational()
  : nm(0), dv(1) { }

inline rational::rational(long n)
  : nm(n), dv(1) { }
inline rational::rational(int n)
  : nm(n), dv(1) { }

inline rational::rational(long n, long d)
  : nm(d < 0 ? -1*n : n), dv(d < 0 ? -1*d : d) { }

 inline rational::rational(double v)
   : nm(0), dv(0)
 {
   *this = dtor(v);
 }

inline rational::rational(const rational& r)
  : nm(r.nm), dv(r.dv) { }

inline rational::rational(const rational::XR& x)
  : nm(x.x_nm), dv(x.x_dv) { };

inline long rational::numerator() const { return nm; }

inline long rational::divisor() const { return dv; }

inline long rational::sign() const
{
  return (nm < 0 ? -1 : (nm > 0 ? 1 : 0));
}

inline bool rational::zero() const
{
  return nm == 0;
}

inline bool rational::finite() const
{
  return dv != 0;
}

inline bool rational::infinite() const
{
  return dv == 0;
}

inline bool rational::integral() const
{
  return dv == 1;
}

inline rational rational::reduce(rational r)
{
  if (r.infinite()) return infinity(r.sign());
  if (r.sign() == 0) return rational(0,1);
  long c = gcd(r.nm, r.dv);
  return rational(r.nm / c, r.dv / c);
}

inline void rational::split(const rational r, long& s, long& k, long& m)
{
  s = r.sign();
  k = ((s * r.nm) / r.dv);
  m = (s * r.nm) - (k * r.dv);
}

inline rational rational::invert(const rational r)
{
  return rational(r.dv, r.nm);
}

inline rational rational::infinity(const long s)
{
  return rational((s < 0 ? -1 : (s > 0 ? 1 : 0)), 0);
}

inline rational rational::infinity(const rational r)
{
  return rational(r.sign(),0);
}

inline rational rational::floor(const rational r)
{
  if (r.infinite()) return r;
  return rational(r.nm / r.dv);
}

inline rational rational::floor_to(const rational r, long d)
{
  if (r.infinite()) return r;
  return rational((r.nm * d) / r.dv, d);
}

inline rational rational::ceil(const rational r)
{
  if (r.infinite()) return r;
  if (r.dv == 1) return r;
  return rational((r.nm / r.dv) + 1, 1);
}

inline rational rational::ceil_to(const rational r, long d)
{
  if (r.infinite()) return r;
  long x = (r.nm * d);
  long y = x / r.dv;
  if ((y * r.dv) == x)
    return rational(y, d);
  else
    return rational(y + 1, d);
}

inline rational rational::frac(const rational r) {
  if (r.infinite()) return r;
  return reduce(rational(r.nm % r.dv, r.dv));
}

inline rational rational::round(const rational r, long d_max)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r << ") round " << d_max << ::std::endl;
#endif
  if (r.infinite()) return r;
  rational s(r);
  while (s.dv > d_max) {
    s.dv = (s.dv / 2);
    s.nm = (s.nm / 2);
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
    ::std::cerr << "s = " << s.nm << " / " << s.dv << ::std::endl;
#endif
    s.reduce();
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
    ::std::cerr << "s = " << s << ::std::endl;
#endif
  }
  return s;
}

inline rational rational::rgcd(const rational r0, const rational r1)
{
  long c = gcd(r0.divisor(), r1.divisor());
  long a0 = r0.numerator() * (r1.divisor() / c);
  long a1 = r1.numerator() * (r0.divisor() / c);
  long d = gcd(a0, a1);
  return rational(d, (r0.divisor() / c) * (r1.divisor() / c) * c).reduce();
}

inline rational rational::min(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    if (r0.sign() < 0) return r0;
    else return r1;
  }
  else if (r1.infinite()) {
    if (r1.sign() < 0) return r1;
    else return r0;
  }
  else if (r0 < r1) return r0;
  else return r1;
}

inline rational rational::max(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    if (r0.sign() > 0) return r0;
    else return r1;
  }
  else if (r1.infinite()) {
    if (r1.sign() > 0) return r1;
    else return r0;
  }
  else if (r0 < r1) return r1;
  else return r0;
}

inline rational rational::reduce() const
{
  return reduce(*this);
}

inline void rational::split(long& s, long& k, long& m) const
{
  split(*this, s, k, m);
}

inline rational rational::invert() const
{
  return invert(*this);
}

inline rational rational::floor() const
{
  return floor(*this);
}

inline rational rational::floor_to(long d) const
{
  return floor_to(*this, d);
}

inline rational rational::ceil() const
{
  return ceil(*this);
}

inline rational rational::frac() const
{
  return frac(*this);
}

inline rational rational::round(long d_max) const
{
  return round(*this, d_max);
}

inline rational rational::round() const
{
  return round(*this, SAFE_RATIONAL_PRECISION);
}

inline rational rational::operator=(const rational r)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational = rational(" << r << ")" << ::std::endl;
#endif
  nm = r.nm;
  dv = r.dv;
  return *this;
}

inline rational rational::operator=(long n)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational = long(" << n << ")" << ::std::endl;
#endif
  nm = n;
  dv = 1;
  return *this;
}

inline rational rational::operator+=(const rational r)
{
  return *this = (*this + r);
}

inline rational rational::operator-=(const rational r)
{
  return *this = (*this - r);
}

inline rational rational::operator*=(const rational r)
{
  return *this = (*this * r);
}

inline rational rational::operator/=(const rational r)
{
  return *this = (*this / r);
}

inline rational rational::operator+=(long n)
{
  return *this = (*this + n);
}

inline rational rational::operator-=(long n)
{
  return *this = (*this - n);
}

inline rational rational::operator*=(long n)
{
  return *this = (*this * n);
}

inline rational rational::operator/=(long n)
{
  return *this = (*this / n);
}

inline double rational::decimal() const { //std::cout<<" n"<<nm<<" "<<dv<<"d "<<std::endl;
return nm/(double)dv; };

inline bool operator==(const rational r0, const rational r1)
{
  if (r0.infinite() && r1.infinite()) return r0.sign() == r1.sign();
  else return ((r0.numerator() == r1.numerator()) &&
	       (r0.divisor() == r1.divisor()));
}

inline bool operator==(const rational r0, long n1)
{
  return ((r0.numerator() == n1) && (r0.divisor() == 1));
}

inline bool operator==(long n0, const rational r1)
{
  return ((r1.numerator() == n0) && (r1.divisor() == 1));
}

inline bool operator!=(const rational r0, const rational r1)
{
  return !(r0 == r1);
}

inline bool operator!=(const rational r0, long n1)
{
  return !(r0 == n1);
}

inline bool operator!=(long n0, const rational r1)
{
  return !(n0 == r1);
}

inline bool operator<(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    // r0 == -INF: r0 < r1 iff r1 != -INF
    if (r0.sign() < 0) {
      if (r1.infinite() && (r1.sign() < 0))
	return false;
      else
	return true;
    }
    // r0 == +INF: r0 < r1 is FALSE
    return false;
  }
  else if (r1.infinite()) {
    // r0 is finite and r1 == -INF: r0 < r1 is FALSE
    if (r1.sign() < 0)
      return false;
    // r0 is finite and r1 == +INF: r0 < r1 is TRUE
    return true;
  }
  // both are finite
  else {
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
    if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
      return (r0.numerator() < r1.numerator());
    }
#endif
    long s0, k0, m0, s1, k1, m1;
    r0.split(s0, k0, m0);
    r1.split(s1, k1, m1);
    // if r0 is negative
    if (s0 < 0) {
      // if r0 is negative and r1 zero or positive, r0 < r1 is TRUE
      if (s1 >= 0)
	return true;
      // both are negative: r0 < r1 iff k0 + (m0/r1.dv) > k1 + (m1/r1.dv)
      else {
	// if r0 (neg) integral part > r1 (neg) integral, r0 < r1 is TRUE
	if (k0 > k1) return true;
	// if r0 (neg) integral part < r1 (neg) integral, r0 < r1 is FALSE
	if (k0 < k1) return false;
	// integral parts and signs are same, have to compare fractions
	return ((rational(m0, r0.divisor()) -
		 rational(m1, r1.divisor())).sign() > 0);
      }
    }
    // r0 is zero or positive
    else {
      // if r1.sign < r0.sign, r0 < r1 is FALSE
      if (s1 < s0)
	return false;
      // both are zero or pos: r0 < r1 iff k0 + (m0/r1.dv) < k1 + (m1/r1.dv)
      else {
	// if r0 integral part < r1 integral, r0 < r1 is TRUE
	if (k0 < k1) return true;
	// if r0 integral part > r1 integral, r0 < r1 is FALSE
	if (k0 > k1) return false;
	// integral parts and signs are same, have to compare fractions
	return ((rational(m0, r0.divisor()) -
		 rational(m1, r1.divisor())).sign() < 0);
      }
    }
  }
}

inline bool operator<=(const rational r0, const rational r1)
{
  if (r0 == r1) return true;
  return (r0 < r1);
}

inline bool operator>(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    // r0 == +INF: r0 > r1 iff r1 != +INF
    if (r0.sign() >= 0) {
      if (r1.infinite() && (r1.sign() >= 0))
	return false;
      else
	return true;
    }
    // r0 == -INF: r0 > r1 is FALSE
    return false;
  }
  else if (r1.infinite()) {
    // r0 is finite and r1 == +INF: r0 > r1 is FALSE
    if (r1.sign() >= 0)
      return false;
    // r0 is finite and r1 == -INF: r0 > r1 is TRUE
    return true;
  }
  // both are finite
  else {
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
    if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
      return (r0.numerator() > r1.numerator());
    }
#endif
    long s0, k0, m0, s1, k1, m1;
    r0.split(s0, k0, m0);
    r1.split(s1, k1, m1);
    // if r0 is negative
    if (s0 < 0) {
      // if r0 is negative and r1 zero or positive, r0 > r1 is FALSE
      if (s1 >= 0)
	return false;
      // both are negative: r0 > r1 iff k0 + (m0/r1.dv) < k1 + (m1/r1.dv)
      else {
	// if r0 (neg) integral part < r1 (neg) integral, r0 > r1 is TRUE
	if (k0 < k1) return true;
	// if r0 (neg) integral part > r1 (neg) integral, r0 > r1 is FALSE
	if (k0 > k1) return false;
	// integral parts and signs are same, have to compare fractions
	return ((rational(m0, r0.divisor()) -
		 rational(m1, r1.divisor())).sign() < 0);
      }
    }
    // r0 is zero or positive
    else {
      // if r1.sign < r0.sign, r0 > r1 is TRUE
      if (s1 < s0)
	return true;
      // both are zero or pos: r0 > r1 iff k0 + (m0/r1.dv) > k1 + (m1/r1.dv)
      else {
	// if r0 integral part > r1 integral, r0 > r1 is TRUE
	if (k0 > k1) return true;
	// if r0 integral part < r1 integral, r0 > r1 is FALSE
	if (k0 < k1) return false;
	// integral parts and signs are same, have to compare fractions
	return ((rational(m0, r0.divisor()) -
		 rational(m1, r1.divisor())).sign() > 0);
      }
    }
  }
}

inline bool operator>=(const rational r0, const rational r1)
{
  if (r0 == r1) return true;
  return (r0 > r1);
}

inline rational operator+(const rational r0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") + rational(" << r1 << ")"
	    << ::std::endl;
#endif
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() == r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " + " << r1 << " not defined"
	          << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() + r1.numerator(), 1);
  }
#endif
  else {
    long c = gcd(r0.divisor(), r1.divisor());
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long long c0 = r0.numerator();
    c0 *= (r1.divisor() / c);
    long long c1 = r1.numerator();
    c1 *= (r0.divisor() / c);
    long long c2 = (c0 + c1);
    long long c3 = r0.divisor();
    c3 *= (r1.divisor() / c);
    if ((c2 > LONG_MAX) || (c2 < LONG_MIN) ||
	(c3 > LONG_MAX) || (c3 < LONG_MIN)) {
      ::std::cerr << "error: numeric overflow in " << r0 << " + " << r1
		  << " (gcd = " << c << ", num = " << c2 << ", div = "
		  << c3 << ")" << ::std::endl;
      abort();
    }
#endif
    long n = (r0.numerator() * (r1.divisor() / c) +
	      r1.numerator() * (r0.divisor() / c));
    long d = ((r0.divisor() / c) * r1.divisor());
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
    ::std::cerr << "c = " << c << ", r0.n = " << r0.numerator()
		<< ", r1.d / c = " << r1.divisor() / c
		<< ", r1.n = " << r1.numerator()
		<< ", r0.d / c = " << r0.divisor() / c
		<< ", result = " << r0.numerator() * (r1.divisor() / c)
		<< " + " << r1.numerator() * (r0.divisor() / c)
		<< " = "  << n << " / " << d
		<< ::std::endl;
#endif
    return rational(n, d).reduce();
  }
}

inline rational operator-(const rational r0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") - rational(" << r1 << ")"
	      << ::std::endl;
#endif
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() != r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " - " << r1 << " not defined"
		    << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1.sign() * -1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() - r1.numerator(), 1);
  }
#endif
  else {
    long c = gcd(r0.divisor(), r1.divisor());
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long long c0 = r0.numerator();
    c0 *= (r1.divisor() / c);
    long long c1 = r1.numerator();
    c1 *= (r0.divisor() / c);
    long long c2 = (c0 - c1);
    long long c3 = r0.divisor();
    c3 *= (r1.divisor() / c);
    if ((c2 > LONG_MAX) || (c2 < LONG_MIN) ||
	(c3 > LONG_MAX) || (c3 < LONG_MIN)) {
      ::std::cerr << "error: numeric overflow in " << r0 << " - " << r1
		  << " (gcd = " << c << ", num = " << c2 << ", div = "
		  << c3 << ")" << ::std::endl;
      abort();
    }
#endif
    long n = (r0.numerator() * (r1.divisor() / c) -
	      r1.numerator() * (r0.divisor() / c));
    long d = ((r0.divisor() / c) * r1.divisor());
    return rational(n, d).reduce();
  }
}

inline rational operator*(const rational r0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") * rational(" << r1 << ")"
	      << ::std::endl;
#endif
  if (r0.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
  else if (r1.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() * r1.numerator(), 1);
  }
#endif
  else {
    long c0 = gcd(r0.numerator(), r1.divisor());
    long c1 = gcd(r1.numerator(), r0.divisor());
    return rational((r0.numerator() / c0) * (r1.numerator() / c1),
		    (r0.divisor() / c1) * (r1.divisor() / c0)).reduce();
  }
}

inline rational operator/(const rational r0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") / rational(" << r1 << ")"
	      << ::std::endl;
#endif
  return (r0 * r1.invert());
}

inline rational operator+(const rational r0, long n1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") + long(" << n1 << ")"
	    << ::std::endl;
#endif
  return rational(r0.numerator() + (n1 * r0.divisor()), r0.divisor()).reduce();
}

inline rational operator-(const rational r0, long n1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") - long(" << n1 << ")"
	    << ::std::endl;
#endif
  return rational(r0.numerator() - (n1 * r0.divisor()), r0.divisor()).reduce();
}

inline rational operator*(const rational r0, long n1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") * long(" << n1 << ")"
	    << ::std::endl;
#endif
  return rational(r0.numerator() * n1, r0.divisor()).reduce();
}

inline rational operator/(const rational r0, long n1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") / long(" << n1 << ")"
	    << ::std::endl;
#endif
  return rational((n1 < 0 ? r0.numerator() * -1 : r0.numerator()),
		  r0.divisor() * (n1 < 0 ? n1 * -1 : n1)).reduce();
}

inline rational operator+(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") + rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return rational(r1.numerator() + (n0 * r1.divisor()), r1.divisor()).reduce();
}

inline rational operator-(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") - rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return rational((n0 * r1.divisor()) - r1.numerator(), r1.divisor()).reduce();
}

inline rational operator*(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") * rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return rational(r1.numerator() * n0, r1.divisor()).reduce();
}

inline rational operator/(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") / rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return (n0 * r1.invert());
}

inline rational safeadd(const rational r0, const rational r1)
{
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() == r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " + " << r1 << " not defined"
		    << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
  else {
    long c = gcd(r0.divisor(), r1.divisor());
    long n0 = r0.numerator();
    long d0 = r0.divisor() / c;
    long n1 = r1.numerator();
    long d1 = r1.divisor() / c;
    assert(d0 > 0);
    assert(d1 > 0);
    while (((LONG_MAX / (2*d1)) < (imag(n0) + 1)) ||
	   ((LONG_MAX / (2*d0)) < (imag(n1) + 1)) ||
	   ((LONG_MAX / (d0 * c)) < (d1 + 1))) {
      if ((d1 < 2) && (d0 < 2)) {
	std::cerr << "error: overflow in safeadd(" << r0 << ", " << r1 << ")"
		  << std::endl;
	abort();
      }
      if (d0 < 2) {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
      else if (d1 < 2) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else if ((imag(n0) - d0) > (imag(n1) - d1)) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
    }
    long f0 = n0 * d1;
    long f1 = n1 * d0;
    long n = f0 + f1;
    long d = d0 * d1 * c;
    return rational(n, d).reduce();
  }
}

inline rational safemul(const rational r0, const rational r1)
{
  if (r0.infinite() || r1.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
  else {
    long c0 = gcd(r0.numerator(), r1.divisor());
    long c1 = gcd(r1.numerator(), r0.divisor());
    long n0 = r0.numerator() / c0;
    long n1 = r1.numerator() / c1;
    long d0 = r0.divisor() / c1;
    long d1 = r1.divisor() / c0;
    // while (((LONG_MAX / imag(n0)) < (imag(n1) + 1)) ||
    //	   ((LONG_MAX / d0) < (d1 + 1)))
    while (((ilog(n0) + ilog(n1)) > LONG_BITS) ||
	   ((ilog(d0) + ilog(d1)) > LONG_BITS)) {
      if ((d1 < 2) && (d0 < 2)) {
	std::cerr << "error: overflow in safemul(" << r0 << ", " << r1 << ")"
		  << std::endl;
	abort();
      }
      if (d0 < 2) {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
      else if (d1 < 2) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else if (ilog(d0) > ilog(d1)) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
    }
    return rational(n0 * n1, d0 * d1).reduce();
  }
}

inline ::std::ostream& operator<<(::std::ostream& s, const rational r)
{
  if (r.infinite()) {
    if (r.sign() < 0) return s << "-INF";
    else return s << "INF";
  }
  else if (r.integral()) {
    return s << r.numerator();
  }
  else {
    return s << r.numerator() << '/' << r.divisor();
  }
}

inline ::std::ostream& operator<<
(::std::ostream& s, const rational::printdecimal& d)
{
  if (d.r.infinite()) {
    if (d.r.sign() < 0) return s << "-INF";
    else return s << "INF";
  }
  else {
    rational ipart = d.r.floor();
    rational dpart = d.r.sign() * d.r.frac();
    s << ipart;
    if (!dpart.zero()) {
      s << '.';
      unsigned int n = 0;
      while (!dpart.zero() && (n < RATIONAL_PRINT_MAX_DECIMALS)) {
	dpart = safemul(dpart, rational(10,1));
	rational digit = dpart.floor();
	s << digit;
	dpart = dpart.frac();
	n += 1;
      }
    }
    return s;
  }
}

END_HSPS_NAMESPACE

#endif
