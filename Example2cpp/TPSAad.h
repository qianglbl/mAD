/*
!****************************
!
!*** Copyright Notice ***
!
!A multi-language auto differentiation package (mAD) Copyright (c) 2025, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
!
!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Intellectual Property Office at
!IPO@lbl.gov.
!
!NOTICE.  This Software was developed under funding from the U.S. Department
!of Energy and the U.S. Government consequently retains certain rights.  As
!such, the U.S. Government has been granted for itself and others acting on
!its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
!Software to reproduce, distribute copies to the public, prepare derivative
!works, and perform publicly and display publicly, and to permit others to do so.
!****************************
! C++ class for forward auto-differentiation
! Ji Qiang - LBNL, jqiang@lbl.gov

*
*/

#ifndef TPS
#define TPS

#include <vector>
#include <iostream>

using namespace std;

class TPSAad{
//private:
//  int terms;
public:
  //dimmax: maximum # of variables to be differentiable
  //For better performance, this number should be the # differentiable variables.
  static const int dimmax = 7;
  double map[dimmax+1] = {0.0,};
  int terms = dimmax+1;

  inline const int getterms() const {return terms;}
  inline const int getnvar() const {return terms-1;}
  void assign(const double a);
  void assign(const double a, const int &ivar);

  TPSAad();
  TPSAad(const int &nvar);
  TPSAad(const double &a);
  TPSAad(const double &a,const int &ivar);
  TPSAad(const TPSAad &a);

  TPSAad& operator=(const double &a);
  TPSAad& operator=(const TPSAad &M);

  TPSAad& operator+=(const TPSAad &);
  TPSAad& operator+=(const double &);
  TPSAad& operator+=(const int &);
  inline const TPSAad operator+(const TPSAad &M) const {return TPSAad(*this)+=M;}
  inline const TPSAad operator+(const double &a) const {return TPSAad(*this)+=a;}
  inline const TPSAad operator+(const int &a) const {return TPSAad(*this)+=a;}
  const friend TPSAad operator+(const int &a,const TPSAad &M);
  const friend TPSAad operator+(const double &a,const TPSAad &M);

  TPSAad& operator-=(const TPSAad &);
  TPSAad& operator-=(const double &);
  TPSAad& operator-=(const int &);
  inline const TPSAad operator-(const TPSAad &M) const {return TPSAad(*this)-=M;}
  inline const TPSAad operator-(const double &a) const {return TPSAad(*this)-=a;}
  inline const TPSAad operator-(const int &a) const {return TPSAad(*this)-=a;}
  const friend TPSAad operator-(const int &a,const TPSAad &M);
  const friend TPSAad operator-(const double &a,const TPSAad &M);

  TPSAad& operator*=(const TPSAad &);
  TPSAad& operator*=(const double &);
  TPSAad& operator*=(const int &);
  inline const TPSAad operator*(const TPSAad &M) const {return TPSAad(*this)*=M;}
  inline const TPSAad operator*(const double &a) const {return TPSAad(*this)*=a;}
  inline const TPSAad operator*(const int &a) const {return TPSAad(*this)*=a;}
  const friend TPSAad operator*(const int &a,const TPSAad &M);
  const friend TPSAad operator*(const double &a,const TPSAad &M);

  TPSAad& operator/=(const TPSAad &);
  TPSAad& operator/=(const double &);
  TPSAad& operator/=(const int &);
  inline const TPSAad operator/(const TPSAad &M) const {return TPSAad(*this)/=M;}
  inline const TPSAad operator/(const double &a) const {return TPSAad(*this)/=a;}
  inline const TPSAad operator/(const int &a) const {return TPSAad(*this)/=a;}
  const friend TPSAad operator/(const int &a,const TPSAad &M);
  const friend TPSAad operator/(const double &a,const TPSAad &M);

  const bool operator==(const TPSAad &M);

  inline const TPSAad operator-() const {return TPSAad(*this)*=(-1);}
  inline const TPSAad operator+() const {return TPSAad(*this);}

  inline const double tpsaval() const {return map[0];}

  friend TPSAad inv(const TPSAad &);
  friend TPSAad exp(const TPSAad &);
  friend TPSAad log(const TPSAad &);
  friend TPSAad sqrt(const TPSAad &);
  friend TPSAad pow(const TPSAad &, const double &);
  friend TPSAad pow(const TPSAad &, const int &);
  friend TPSAad sin(const TPSAad &);
  friend TPSAad cos(const TPSAad &);
  friend TPSAad tan(const TPSAad &);
  friend TPSAad asin(const TPSAad &);
  friend TPSAad acos(const TPSAad &);
  friend TPSAad atan(const TPSAad &);
  friend TPSAad sinh(const TPSAad &);
  friend TPSAad cosh(const TPSAad &);
  friend TPSAad tanh(const TPSAad &);
  friend TPSAad asinh(const TPSAad &);
  friend TPSAad acosh(const TPSAad &);
  friend TPSAad atanh(const TPSAad &);

};



#endif

