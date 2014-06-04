#include <math.h>
#include <iostream>
#include "Vector.hpp"
#define ROUND_CUTOFF 0.00001 

//This file uses some conventions that do NOT optimize for speed; rather, for
//flexibility. This should be relatively easy to adapt to 3d, for example. 

//Vector::Dimension = dim;

void roundFloat(float_T& f, const float_T precision=ROUND_CUTOFF) {
  if (fabs(round(f)-f) < precision) {
    f = round(f);
    if (f == -0)
      f = 0;
  }
}

Vector::Vector() {
  for (size_t i = 0; i < Vector::Dimension; i++)
    data_[i] = 0;
}

Vector::Vector(const Vector& copy) 
{
  for (size_t i = 0; i < Vector::Dimension; i++)
    data_[i] = copy[i];
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

bool Vector::isZero() const {
  return (this->norm() == 0);
}

Vector Vector::projectOn(const Vector& base) const {
  if (base.isZero())
    return base;
  else 
    return ((eucInnerProd(base, (*this))/eucInnerProd(base,base)) * base);
}

void Vector::normalize() {
  float_T mag = this->norm();
  for (size_t i=0; i < Vector::Dimension; ++i) 
    (*this)[i] *= 1/mag;
}
