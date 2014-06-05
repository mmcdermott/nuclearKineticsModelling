#include <math.h>
#include <iostream>
#include <assert.h>
#include "Vector.hpp"
#define ROUND_CUTOFF 0.00001 

//This file uses some conventions that do NOT optimize for speed; rather, for
//flexibility. This should be relatively easy to adapt to 3d, for example. 

void roundFloat(float_T& f, const float_T precision=ROUND_CUTOFF) {
  /* roundFloat: A function to abstract away rounding a float to an arbitrary
   * precision.
   * Inputs:
   *   f: The float to be rounded.
   *   precision: the precision to which f should be rounded.
   * Output: f rounded to precision.
   */
  if (fabs(round(f)-f) < precision) {
    f = round(f);
    if (f == -0)
      f = 0;
  }
}

//Vector Class Functions: 
// Constructors:
Vector::Vector() {
  /* Vector::Vector(): A default constructor. 
   * Inputs: (none)
   * Output: The zero vector.
   */
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = 0;
}

Vector::Vector(const Vector& copy) {
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = copy[i];
}

Vector::Vector(const std::vector<float_T>& copy) {
  assert(copy.size() == Vector::Dimension);
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = copy[i];
}

Vector& Vector::operator=(const Vector& copy) {
  for (size_t i = 0; i < Vector::Dimension; i++)
    (*this)[i] = copy[i];
  return (*this);
}

const float_T& Vector::operator[](const unsigned int i) const {
  return this->data_[i];
}

float_T& Vector::operator[](const unsigned int i) {
  return this->data_[i];
}

float_T Vector::norm() const {
  float_T result = 0;
  for (size_t i=0; i < Vector::Dimension; ++i)
    result += pow((*this)[i],2);
  result = sqrt(result);
  roundFloat(result);
  return result;
}

void Vector::round() {
  for (size_t i = 0; i < Vector::Dimension; ++i) {
    roundFloat((*this)[i]);
  }
}

float_T eucInnerProd(const Vector& v1, const Vector& v2) {
  float_T result = 0;
  for (size_t i=0; i < Vector::Dimension; ++i) {
    result += v1[i] * v2[i];
  }
  roundFloat(result);
  return result;
}

bool orthogonal(const Vector& v1, const Vector& v2) {
  return (eucInnerProd(v1, v2) == 0);
}

Vector operator*(const float_T s, const Vector& v) {
  Vector result;
  for (size_t i = 0; i < Vector::Dimension; ++i)
    result[i] = s*v[i];
  return result;
}

Vector Vector::operator+(const Vector& rhs) const {
  Vector result;
  for (size_t i = 0; i < Vector::Dimension; ++i)
    result[i] = (*this)[i] + rhs[i];
  return result;
}

Vector Vector::operator-(const Vector& rhs) const {
  return (*this) + ((-1)*rhs);
}

void Vector::operator+=(const Vector& rhs) {
  (*this) = (*this) + rhs;
}

void Vector::operator-=(const Vector& rhs) {
  (*this) = (*this) - rhs;
}

void Vector::operator *=(const float_T s) {
  (*this) = s*(*this);
}

bool Vector::isZero() const {
  return (this->norm() == 0);
}

void Vector::zero() {
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = 0;
}


Vector Vector::projectOn(const Vector& base) const {
  if (base.isZero())
    return base;
  else 
    return ((eucInnerProd(base, (*this))/eucInnerProd(base,base)) * base);
}

Vector Vector::normalize() const {
  Vector result = (*this);
  float_T mag = result.norm();
  result *= (1/mag);
  return result;
}

float_T distBetween(const Vector& a, const Vector& b) {
  return (a-b).norm();
}

Vector operator-(const Vector& v) {
  return (-1)*v;
}
