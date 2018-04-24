class Complex{
private:
  double real;
  double imag;
public:
  Complex(double real=0.0, double imag=0.0)
            :m_real(real), m_imag(imag) { };/////////
  Complex Add(Complex A, Complex B);

};
Complex Complex::Add(Complex A, Complex B){
  Complex res;
  res.real = A.real+B.real;
  res.imag = A.imag+B.imag;
  return res;
}

class FFT_v1{
private:

public:  

};


/* reference
#include <iostream>
using namespace std;

template<typename T>
class Complex {
    private:
        T m_real;
        T m_imag;
    public:
        Complex<T>(T real=0.0, T imag=0.0)
            :m_real(real), m_imag(imag) { };
        Complex<T>(const Complex<T>& c)
            :m_real(c.m_real), m_imag(c.m_imag) { };
        Complex<T> operator -() const;
        Complex<T>& operator +=(const Complex<T>&);
        Complex<T>& operator -=(const Complex<T>&);
        Complex<T>& operator *=(const Complex<T>&);
        Complex<T>& operator /=(const Complex<T>&);
        Complex<T> operator +(const Complex<T>&) const;
        Complex<T> operator -(const Complex<T>&) const;
        Complex<T> operator *(const Complex<T>&) const;
        Complex<T> operator /(const Complex<T>&) const;
        template<typename U>
        friend ostream& operator <<(ostream&, const Complex<U>&);
};

template<typename T>
Complex<T> Complex<T>::operator -() const {
    return Complex(-m_real, -m_imag);
}

template<typename T>
Complex<T>& Complex<T>::operator +=(const Complex<T>& c2) {
    m_real += c2.m_real;
    m_imag += c2.m_imag;
    return *this;
}

template<typename T>
Complex<T>& Complex<T>::operator -=(const Complex<T>& c2) {
    return *this += -c2;
}

template<typename T>
Complex<T>& Complex<T>::operator *=(const Complex<T>& c2) {
    T real = m_real * c2.m_real - m_imag * c2.m_imag;
    T imag = m_real * c2.m_imag + m_imag * c2.m_real;
    m_real = real;
    m_imag = imag;
    return *this;
}

template<typename T>
Complex<T>& Complex<T>::operator /=(const Complex<T>& c2) {
    Complex<T> nm = Complex<T>(*this) * Complex<T>(c2.m_real, -c2.m_imag);
    T dn = c2.m_real * c2.m_real + c2.m_imag * c2.m_imag;
    m_real = nm.m_real / dn;
    m_imag = nm.m_imag / dn;
    return *this;
}

template<typename T>
Complex<T> Complex<T>::operator +(const Complex<T>& c2) const {
    return Complex<T>(*this) += c2;
}

template<typename T>
Complex<T> Complex<T>::operator -(const Complex<T>& c2) const {
    return Complex<T>(*this) -= c2;
}

template<typename T>
Complex<T> Complex<T>::operator *(const Complex<T>& c2) const {
    return Complex<T>(*this) *= c2;
}

template<typename T>
Complex<T> Complex<T>::operator /(const Complex<T>& c2) const {
    return Complex<T>(*this) /= c2;
}

template<typename T>
ostream& operator <<(ostream& os, const Complex<T>& c) {
    os << c.m_real;
    if (c.m_imag < 0)
        os << c.m_imag << "i";
    else if (c.m_imag > 0)
        os << "+" << c.m_imag << "i";
    return os;
}
*/
