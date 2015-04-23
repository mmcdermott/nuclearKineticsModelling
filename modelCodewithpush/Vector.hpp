#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED 1

#include <math.h>
#include <vector>

//Abstracting the level of precsision used:
typedef double float_T;
const int dim = 2;

class Vector {
  /* The Vector class enables 2d and 3d vector arithmetic, for abstraction
   * purposes within the MT modeling code. 
   */
  public: 
    //Constructors: 
    Vector();
    Vector(const Vector& copy);
    Vector(const std::vector<float_T>& copy);
    
    //Assignment operator:
    Vector& operator=(const Vector& copy);
    //Access operators (const and non-const versions):
    const float_T& operator[](const unsigned int i) const;
    float_T& operator[](const unsigned int i);

    //Arithmetic operators:
    Vector operator+(const Vector& rhs) const;
    Vector operator-(const Vector& rhs) const;
    void operator+=(const Vector& rhs);
    void operator-=(const Vector& rhs);
    void operator*=(const float_T s);

    //Vector operators:
    float_T norm() const;
    Vector projectOn(const Vector& base) const;
    Vector normalize() const;
    //Misc. Operators:
    void round();
    bool isZero() const;
    void zero();

    //The Dimension of the vectors being used.
    const static unsigned int Dimension = dim;
  private:
    //The data:
    float_T data_[Vector::Dimension];   
    //float_T data_[dim];   
};

//Scalar multicplication and negation:
Vector operator-(const Vector& v);
Vector operator*(const float_T s, const Vector& v);
//Multi-vector functions:
float_T eucInnerProd(const Vector& v1, const Vector& v2);
bool orthogonal(const Vector& v1, const Vector& v2);
float_T distBetween(const Vector& a, const Vector& b);

#endif
