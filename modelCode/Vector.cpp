#include <math.h>
#include <iostream>
#include <assert.h>
#include "Vector.hpp"
#define ROUND_CUTOFF 0.00001 

//This file uses some conventions that do NOT optimize for speed; rather, for
//flexibility. This should be very easy to adapt to 3D, for example (which might
//just take updating Vector::Dimension).

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
  /* A default constructor. 
   * Inputs: (none)
   * Output: The zero vector.
   */
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = 0;
}

Vector::Vector(const Vector& copy) {
  /* A copy constructor.
   * Inputs: 
   *   const Vector& copy: The vector to copy.
   * Output: A duplicate of the input vector, copy. 
   */
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = copy[i];
}

Vector::Vector(const std::vector<float_T>& copy) {
  /* A vector<float_T> constructor.
   * Inputs: 
   *   const std::vector<float_T>& copy): This is the standard vector class
   *     object we want to copy into our class. 
   * Outputs: A duplicate of the input vector, copy.
   */

  //We want to make sure we're not trying to copy a 3D vector to a 2D or
  //something. Occasionally, this will make things break if one tries to do
  //something funny, like initialize a Vector out of a really long std::vector,
  //for some good reason, but in that case it is better stylistically to isolate
  //the components you want a 2 component vector first and build the Vector that
  //way. 
  assert(copy.size() == Vector::Dimension);
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = copy[i];
}

Vector& Vector::operator=(const Vector& copy) {
  /* An assignment operator.
   * Inputs: 
   *   const Vector& copy: The vector to copy. 
   * Output: This technically outputs a reference to a vector, which is really
   *   just a programmatic necessity to make it update the object in question,
   *   such that if we say Vector v = Vector w, then v will become a duplicate
   *   of w. Technically, a reference to v, which has been updated to look like
   *   w, will be output. 
   */
  for (size_t i = 0; i < Vector::Dimension; i++)
    (*this)[i] = copy[i];
  return (*this);
}

const float_T& Vector::operator[](const unsigned int i) const {
  /* A read-only access operator. This enables reading components of a const
   *   Vector. 
   * Inputs: 
   *   const unsigned int i: This is the index of the vector that we want to
   *     access. 
   * Output: The ith component of the vector. 
   */

  //i better be a valid index, else we'll get a segfault. 
  assert(i < Vector::Dimension);
  return this->data_[i];
}

float_T& Vector::operator[](const unsigned int i) {
  /* A read-write access operator. This enables reading and writing components
   *   of a Vector as accessed via their index. 
   * Inputs: 
   *   const unsigned int i: This is the index of the vector that we want to
   *     access. 
   * Output: A modifiable reference to the ith component of the vector. 
   */

  //i better be a valid index, else we'll get a segfault. 
  assert(i < Vector::Dimension);
  return this->data_[i];
}

float_T Vector::norm() const {
  /* A function that returns the euclidean norm of calling Vector.  
   * Inputs: (none)
   * Output: The euclidean norm of the calling Vector.
   */
  float_T result = 0;
  for (size_t i=0; i < Vector::Dimension; ++i)
    result += pow((*this)[i],2);
  result = sqrt(result);
  // We always round our floats here to a particular precision to avoid
  // accumulation of round-off error (which isn't a problem right now, but this
  // should prevent it from becoming a problem).
  roundFloat(result);
  return result;
}

void Vector::round() {
  /* A function that rounds the calling object to the hard-coded
   *   ROUND_PRECISION. I don't believe this is ever used right now. 
   * Inputs: (none)
   * Output: (none, but updates the calling object)
   */
  for (size_t i = 0; i < Vector::Dimension; ++i) {
    roundFloat((*this)[i]);
  }
}

float_T eucInnerProd(const Vector& v1, const Vector& v2) {
  /* A function that computes the euclidean inner product of the input Vectors.
   * Inputs: 
   *   const Vector& v1: The first input vector which will be used to compute
   *     the inner product. 
   *   const Vector& v2: The second input vector which will be used to compute
   *     the inner product. 
   * Output: The euclidean inner product of v1 and v2. 
   */
  float_T result = 0;
  for (size_t i=0; i < Vector::Dimension; ++i) {
    result += v1[i] * v2[i];
  }
  roundFloat(result);
  return result;
}

bool orthogonal(const Vector& v1, const Vector& v2) {
  /* A function to test if the input Vectors are orthogonal. 
   * Inputs: 
   *   const Vector& v1: One of the two input vectors to be tested for
   *     orthogonality. 
   *   const Vector& v2: One of the two input vectors to be tested for
   *     orthogonality. 
   * Output: A boolean answering the question: 'Are v1 and v2 orthogonal?'.
   */
  return (eucInnerProd(v1, v2) == 0);
}

Vector operator*(const float_T s, const Vector& v) {
  /* A function to multiply a Vector by a scalar. 
   * Inputs: 
   *   const float_T s: The input scalar. 
   *   const Vector& v: The input Vector.
   * Output: The scaled vector s*v;
   */
  Vector result;
  for (size_t i = 0; i < Vector::Dimension; ++i)
    result[i] = s*v[i];
  return result;
}

Vector Vector::operator+(const Vector& rhs) const {
  /* A function to add two vectors.
   * Inputs: 
   *   const Vector& rhs: The Vector to add to the calling object. 
   * Output: The sum of the calling object and rhs. 
   */
  Vector result;
  for (size_t i = 0; i < Vector::Dimension; ++i)
    result[i] = (*this)[i] + rhs[i];
  return result;
}

Vector Vector::operator-(const Vector& rhs) const {
  /* A function to take the difference of two vectors. 
   * Inputs: 
   *   const Vector& rhs: The Vector to subtract from the calling object. 
   * Output: The difference of the calling object and rhs. 
   */
  return (*this) + ((-1)*rhs);
}

void Vector::operator+=(const Vector& rhs) {
  /* A function to perform an in place addition and update of the calling
   *   vector.
   * Inputs: 
   *   const Vector& rhs: The vector to add to the calling object. 
   * Output: (none, but the calling object is updated to be equal to the sum of
   *          itself and rhs)
   */
  (*this) = (*this) + rhs;
}

void Vector::operator-=(const Vector& rhs) {
  /* A function to perform an in place subtraction and update of the calling
   *   vector.
   * Inputs: 
   *   const Vector& rhs: The vector to subtract from the calling object. 
   * Output: (none, but the calling object is updated to be equal to the
   *          difference of itself and rhs)
   */
  (*this) = (*this) - rhs;
}

void Vector::operator *=(const float_T s) {
  /* A function to perform an in place scalar multiplication and update of the
   * calling vector.
   * Inputs: 
   *   const float_T s: The scalar to scale the calling object. 
   * Output: (none, but the calling object is updated to be equal to itself,
   *          scaled by s)
   */
  (*this) = s*(*this);
}

bool Vector::isZero() const {
  /* A function to test if the calling Vector is the zero Vector. 
   * Inputs: (none)
   * Output: 'true' if the calling Vector is zero, and 'false' otherwise. 
   */
  return (this->norm() == 0);
}

void Vector::zero() {
  /* A function to zero out the calling Vector. 
   * Inputs: (none)
   * Output: (none, but the calling Vector is set to zero)
   */
  for (size_t i = 0; i < Vector::Dimension; i++)
    this->data_[i] = 0;
}


Vector Vector::projectOn(const Vector& base) const {
  /* A function to project the calling Vector onto the input Vector. 
   * Inputs: 
   *   const Vector& base: The base vector onto which the calling Vector should
   *     be projected. 
   * Output: The projection of the calling Vector onto base. If base == 0, then
   *   the result is simply the zero vector. 
   */
  if (base.isZero())
    return base;
  else 
    return ((eucInnerProd(base, (*this))/eucInnerProd(base,base)) * base);
}

Vector Vector::normalize() const {
  /* A function that normalizes the calling Vector. 
   * Inputs: (none)
   * Output: (none, but the calling Vector is normalized, or scaled to have unit
   *          norm)
   */
  Vector result = (*this);
  float_T mag = result.norm();
  result *= (1/mag);
  return result;
}

float_T distBetween(const Vector& a, const Vector& b) {
  /* A function that calculates the distance between two Vectors. 
   * Inputs: 
   *   const Vector& a: One of the two input Vectors that determine the
   *     calculated distance. 
   *   const Vector& b: One of the two input Vectors that determine the
   *     calculated distance. 
   * Output: The distance between a and b. 
   */
  return (a-b).norm();
}

Vector operator-(const Vector& v) {
  /* The unitary negation operator for an input Vector. 
   * Inputs: 
   *   const Vector& v: The vector to negate. 
   * Output: The negated Vector -v. 
   */
  return (-1)*v;
}
