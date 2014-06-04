#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED 1

#include <math.h>

typedef double float_T;
const int dim = 2;

class Vector {
  public: 
    Vector();
    Vector(const Vector& copy);
    
    Vector& operator=(const Vector& copy);
    const float_T& operator[](const unsigned int i) const;
    float_T& operator[](const unsigned int i);

    float_T norm() const;
    void round();
    bool isZero() const;
    Vector operator+(const Vector& rhs) const;
    Vector operator-(const Vector& rhs) const;
    void operator+=(const Vector& rhs);
    void operator-=(const Vector& rhs);
    Vector projectOn(const Vector& base) const;
    void normalize();

    const static unsigned int Dimension = dim;
  private:
    float_T data_[Vector::Dimension];   
    //float_T data_[dim];   
};

float_T eucInnerProd(const Vector& v1, const Vector& v2);
bool orthogonal(const Vector& v1, const Vector& v2);
Vector operator*(const float_T s, const Vector& v);

#endif
