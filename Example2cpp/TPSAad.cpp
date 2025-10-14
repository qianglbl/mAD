/*
 * TPSAad.cpp: by Ji Qiang
 *
 */
#include "TPSAad.h"
#include <iostream>
#include <cmath>

using namespace std;

void TPSAad::assign(const double a) {
	terms = dimmax + 1;
	map[0] = a;
	for (int i = 1; i < terms; i++) {map[i] = 0.0;}
}

void TPSAad::assign(const double a, const int& ivar) {
        terms = dimmax+1;
	for (int i = 1; i < terms; i++) {map[i] = 0.0;}
	map[ivar] = 1.0;
	map[0] = a;
}

TPSAad::TPSAad() {
	terms = dimmax+1;
	map[terms] = {0.0, };
}

TPSAad::TPSAad(const int& nvar) {
	terms = nvar+1;
	map[terms] = {0.0, };
}

TPSAad::TPSAad(const double& a) {
	terms = dimmax+1;
        map[terms] = {0.0, };
        map[0] = a;
}

TPSAad::TPSAad(const double& a,const int& ivar) {
	terms = dimmax+1;
        map[terms] = {0.0, };
        map[0] = a;
        map[ivar] = 1.0;
}

TPSAad::TPSAad(const TPSAad &M) {
        terms = M.terms;
        for (int i = 0;i<M.terms;i++) {map[i] = M.map[i];}
}

TPSAad& TPSAad::operator=(const double &a ){
        this->terms=dimmax+1;
        this->map[0]=a;
        return *this;
}

TPSAad& TPSAad::operator=(const TPSAad &M){
        if (this != &M){
                this->terms=M.terms;
	        for (int i = 0; i<M.terms;i++) {this->map[i] = M.map[i];}
        }
        return *this;
}

TPSAad& TPSAad::operator+=(const TPSAad &M) {
	for (int i = 0; i<M.terms; i++) {this->map[i] += M.map[i];}
	return *this;
}

TPSAad& TPSAad::operator+=(const double &a) {
	this->map[0] += a; return *this;
}

TPSAad& TPSAad::operator+=(const int &a) {
	this->map[0] += a; return *this;
}

TPSAad& TPSAad::operator-=(const TPSAad &M) {
	for (int i = 0; i<M.terms; i++) {this->map[i] -= M.map[i];}
	return *this;
}

TPSAad& TPSAad::operator-=(const double &a) {
	this->map[0] -= a; return *this;
}

TPSAad& TPSAad::operator-=(const int &a) {
	this->map[0] -= a; return *this;
}

TPSAad& TPSAad::operator*=(const TPSAad &M) {
	//int multmax = max(M.terms, this->terms);
	for (int i = 1; i<this->terms; i++) {
		this->map[i] = M.map[0]*this->map[i] + this->map[0]*M.map[i];

	} 
	this->map[0] *= M.map[0];
	return *this;
}

TPSAad& TPSAad::operator*=(const double &a) {
	for (int i = 0; i<this->terms; i++) {this->map[i] *= a;} return *this;
}

TPSAad& TPSAad::operator*=(const int &a) {
	for (int i = 0; i<this->terms; i++) {this->map[i] *= a;} return *this;
}


TPSAad& TPSAad::operator/=(const TPSAad &M) {
	if (abs(M.map[0]) < 1.0e-14) { 
		cout << "Error: Divide by zero, in TPSAad" << endl; exit(0);
	}
	for (int i = 1; i<this->terms; i++) {
		this->map[i] = (this->map[i]*M.map[0]-M.map[i]*this->map[0])/(M.map[0]*M.map[0]);
	}
	this->map[0] /= M.map[0];
	return *this;
}

TPSAad& TPSAad::operator/=(const double &a) {
	for (int i = 0; i<this->terms; i++) {
	  this->map[i] /= a; 
	}
	return *this;
}

TPSAad& TPSAad::operator/=(const int &a) {
	for (int i = 0; i<this->terms; i++) {
	  this->map[i] /= a; 
	}
	return *this;
}

  const TPSAad operator+(const int &a, const TPSAad &M)
  {
    TPSAad N=M;
    N.map[0] = M.map[0]+a;
    return N;
  }

  const TPSAad operator+(const double &a, const TPSAad &M)
  {
    TPSAad N=M;
    N.map[0] = M.map[0]+a;
    return N;
  }

  const TPSAad operator-(const int &a, const TPSAad &M)
  {
    TPSAad N=M;
    for (int i = 1; i<M.terms;i++) {N.map[i] = -M.map[i];}
    N.map[0] = a-M.map[0];
    return N;
  }

  const TPSAad operator-(const double &a, const TPSAad &M)
  {
    TPSAad N=M;
    for (int i = 1; i<M.terms;i++) {N.map[i] = -M.map[i];}
    N.map[0] = a-M.map[0];
    return N;
  }

  const TPSAad operator*(const int &a, const TPSAad &M)
  {
    TPSAad N=M;
    for (int i = 0; i<M.terms;i++) {N.map[i] = M.map[i]*a;}
    return N;
  }

  const TPSAad operator*(const double &a, const TPSAad &M)
  {
    TPSAad N=M;
    for (int i = 0; i<M.terms;i++) {N.map[i] = M.map[i]*a;}
    return N;
  }

  const TPSAad operator/(const int &a, const TPSAad &M)
  {
    TPSAad N=M;
    N.map[0] = a/M.map[0];
    for (int i = 1; i<M.terms;i++) { N.map[i] = -a*M.map[i]/(M.map[0]*M.map[0]);}
    return N;
  }

  const TPSAad operator/(const double &a, const TPSAad &M)
  {
    TPSAad N=M;
    N.map[0] = a/M.map[0];
    for (int i = 1; i<M.terms;i++) { N.map[i] = -a*M.map[i]/(M.map[0]*M.map[0]);}
    return N;
  }

const bool TPSAad::operator==(const TPSAad & M) {
	if (this->terms != M.terms || this->map[0] != M.map[0]) {return false;}
	if (sizeof(this->map) != sizeof(M.map)) {return false;}
	if (this->map != M.map) {return false;}
	return true; 
}

TPSAad inv(const TPSAad &M) {
	if (abs(M.map[0]) < 1.0e-14) { 
		cout << "Error: Divide by zero, in TPSAad" << endl; exit(0);
	}
	TPSAad N = M;
	N.map[0] = 1.0/M.map[0];
	for (int i = 1; i<M.terms;i++) { N.map[i] = -M.map[i]/(M.map[0]*M.map[0]);}
	return N;
}

TPSAad exp(const TPSAad &M) {
	TPSAad N = M;
	//TPSAad N;
	double temp = M.map[0];
	N.map[0] = exp(temp);
	for (int i = 1; i<M.terms;i++) {N.map[i] = exp(M.map[0])*M.map[i];}
	return N;
}

TPSAad log(const TPSAad &M) {
	if (abs(M.map[0]) < 1.0e-15) {
		cout << "Zero in log function: " << endl;
		exit(0);
	}
	TPSAad N = M;
	N.map[0] = log(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/M.map[0];}
	return N;
}
TPSAad sqrt(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = sqrt(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = 0.5*M.map[i]/sqrt(M.map[0]);}
	return N;
}

TPSAad pow(const TPSAad &M, const int &a) {
        TPSAad N = M;
        if (a == 1) {
                return N;
        }
        else if (a == 0) {
                return 1.0;
        }
        else {
                N.map[0] = pow(M.map[0],a);
                for (int i = 1; i<M.terms;i++) {N.map[i] = a*pow(M.map[0], (a-1))*M.map[i];}
                return N;
        }
}

TPSAad pow(const TPSAad &M, const double &a) {
	TPSAad N = M;
	if (abs(a-1) < 1e-14) {
		return N;
	}
	else if (abs(a) < 1e-14) {
		return 1.0;
	}
	else {
		N.map[0] = pow(M.map[0],a);
		for (int i = 1; i<M.terms;i++) {N.map[i] = a*pow(M.map[0], (a-1))*M.map[i];}
		return N;
	}
}

TPSAad sin(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = sin(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = cos(M.map[0])*M.map[i];}
	return N;
}

TPSAad cos(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = cos(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = -sin(M.map[0])*M.map[i];}
	return N;
}

TPSAad tan(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = tan(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = (1.0+pow(tan(M.map[0]),2))*M.map[i];}
	return N;
}
TPSAad asin(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = asin(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/sqrt(1-pow(M.map[0],2));}
        return N;
}
TPSAad acos(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = acos(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = -M.map[i]/sqrt(1-pow(M.map[0],2));}
        return N;
}
TPSAad atan(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = atan(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/(1+pow(M.map[0],2));}
        return N;
}
TPSAad sinh(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = sinh(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = cosh(M.map[0])*M.map[i];}
	return N;
}
TPSAad cosh(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = cosh(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = sinh(M.map[0])*M.map[i];}
	return N;
}
TPSAad tanh(const TPSAad &M) {
	TPSAad N = M;
	N.map[0] = tanh(M.map[0]);
	for (int i = 1; i<M.terms;i++) {N.map[i] = (1.0 + pow(tanh(M.map[0]),2))*M.map[i];}
	return N;
}
TPSAad asinh(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = asinh(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/sqrt(1+pow(M.map[0],2));}
        return N;
}
TPSAad acosh(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = acosh(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/sqrt(pow(M.map[0],2)-1);}
        return N;
}
TPSAad atanh(const TPSAad &M) {
        TPSAad N = M;
        N.map[0] = atanh(M.map[0]);
        for (int i = 1; i<M.terms;i++) {N.map[i] = M.map[i]/(1-pow(M.map[0],2));}
        return N;
}
