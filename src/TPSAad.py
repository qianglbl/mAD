#****************************
#
#*** Copyright Notice ***
#
#A multi-language auto differentiation package (mAD) Copyright (c) 2025.
#The Regents of the University of California, through Lawrence Berkeley National
#Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
#All rights reserved.
#
#If you have questions about your rights to use or distribute this software,
#please contact Berkeley Lab's Intellectual Property Office at
#IPO@lbl.gov.
#
#NOTICE.  This Software was developed under funding from the U.S. Department
#of Energy and the U.S. Government consequently retains certain rights.  As
#such, the U.S. Government has been granted for itself and others acting on
#its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
#Software to reproduce, distribute copies to the public, prepare derivative
#works, and perform publicly and display publicly, and to permit others to do so.
#
#****************************
# Python class for forward auto-differentiation
# Ji Qiang with help of Antropic Claude Sonnet - LBNL, jqiang@lbl.gov
#
#****************************

import numpy as np
from math import *

class TPSAad:
    #dimmax: maximum # of variables to be differentiable
    #For better performance, this number should be the # differentiable variables.
    dimmax = 5

    def __init__(self, a=None, ivar=None):
        self.terms = self.dimmax + 1
        self.map = np.zeros(self.terms)
        
        if isinstance(a, (int, float)) and ivar is None:
            self.map[0] = a
        elif isinstance(a, (int, float)) and isinstance(ivar, int):
            self.map[0] = a
            self.map[ivar] = 1.0

    def assign(self, a, ivar=None):
        if ivar is None:
            self.map[0] = a
            self.map[1:] = 0.0
        else:
            self.map[:] = 0.0
            self.map[0] = a
            self.map[ivar] = 1.0

    def __add__(self, other):
        result = TPSAad()
        if isinstance(other, (int, float)):
            result.map = self.map.copy()
            result.map[0] += other
        else:
            result.map = self.map + other.map
        return result

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        result = TPSAad()
        if isinstance(other, (int, float)):
            result.map = self.map.copy()
            result.map[0] -= other
        else:
            result.map = self.map - other.map
        return result

    def __rsub__(self, other):
        result = TPSAad()
        if isinstance(other, (int, float)):
            result.map = -self.map
            result.map[0] = other - self.map[0]
        return result

    def __mul__(self, other):
        result = TPSAad()
        if isinstance(other, (int, float)):
            result.map = self.map * other
        else:
            result.map[0] = self.map[0] * other.map[0]
            for i in range(1, self.terms):
                result.map[i] = other.map[0]*self.map[i] + self.map[0]*other.map[i]
        return result

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        result = TPSAad()
        if isinstance(other, (int, float)):
            if abs(other) < 1e-14:
                raise ValueError("Division by zero")
            result.map = self.map / other
        else:
            if abs(other.map[0]) < 1e-14:
                raise ValueError("Division by zero")
            result.map[0] = self.map[0] / other.map[0]
            for i in range(1, self.terms):
                result.map[i] = (self.map[i]*other.map[0] - other.map[i]*self.map[0])/(other.map[0]**2)
        return result

    def __rtruediv__(self, other):
        result = TPSAad()
        if abs(self.map[0]) < 1e-14:
            raise ValueError("Division by zero")
        result.map[0] = other / self.map[0]
        for i in range(1, self.terms):
            result.map[i] = -other * self.map[i]/(self.map[0]**2)
        return result

    def __neg__(self):
        result = TPSAad()
        result.map = -self.map
        return result

    @property
    def tpsaval(self):
        return self.map[0]

# Mathematical functions
def sin(M):
    result = TPSAad()
    result.map[0] = np.sin(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = np.cos(M.map[0]) * M.map[i]
    return result

def cos(M):
    result = TPSAad()
    result.map[0] = np.cos(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = -np.sin(M.map[0]) * M.map[i]
    return result

def tan(M):
    """
    Tangent of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.tan(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = (1.0 + np.tan(M.map[0])**2) * M.map[i]
    return result

def asin(M):
    """
    Arcsine (inverse sine) of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arcsin(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i] / np.sqrt(1 - M.map[0]**2)
    return result

def acos(M):
    """
    Arccosine (inverse cosine) of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arccos(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = -M.map[i] / np.sqrt(1 - M.map[0]**2)
    return result

def atan(M):
    """
    Arctangent (inverse tangent) of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arctan(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i] / (1 + M.map[0]**2)
    return result

def exp(M):
    result = TPSAad()
    result.map[0] = np.exp(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = np.exp(M.map[0]) * M.map[i]
    return result

def log(M):
    """
    Natural logarithm of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    Raises:
        ValueError: if argument is zero or negative
    """
    if abs(M.map[0]) < 1.0e-15:
        raise ValueError("Zero in log function")

    result = TPSAad()
    result.map = M.map.copy()  # Create a copy of the map
    result.map[0] = np.log(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i]/M.map[0]
    return result

def sqrt(M):
    """
    Square root of a TPSAad object.

    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.sqrt(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = 0.5 * M.map[i] / np.sqrt(M.map[0])
    return result

def pow(M, a):
    """
    Power function for TPSAad object.

    Args:
        M: TPSAad object
        a: integer or float exponent
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()

    if isinstance(a, int):
        if a == 1:
            return result
        elif a == 0:
            result.map[:] = 0.0
            result.map[0] = 1.0
            return result
        else:
            result.map[0] = np.power(M.map[0], a)
            for i in range(1, M.terms):
                result.map[i] = a * np.power(M.map[0], (a-1)) * M.map[i]
            return result

    elif isinstance(a, float):
        if abs(a-1) < 1e-14:
            return result
        elif abs(a) < 1e-14:
            result.map[:] = 0.0
            result.map[0] = 1.0
            return result
        else:
            result.map[0] = np.power(M.map[0], a)
            for i in range(1, M.terms):
                result.map[i] = a * np.power(M.map[0], (a-1)) * M.map[i]
            return result
    else:
        raise TypeError("Exponent must be int or float")

def sinh(M):
    """
    Hyperbolic sine of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.sinh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = np.cosh(M.map[0]) * M.map[i]
    return result

def cosh(M):
    """
    Hyperbolic cosine of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.cosh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = np.sinh(M.map[0]) * M.map[i]
    return result

def tanh(M):
    """
    Hyperbolic tangent of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.tanh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = (1.0 + np.tanh(M.map[0])**2) * M.map[i]
    return result

def asinh(M):
    """
    Inverse hyperbolic sine of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    """
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arcsinh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i] / np.sqrt(1 + M.map[0]**2)
    return result

def acosh(M):
    """
    Inverse hyperbolic cosine of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    Raises:
        ValueError: if argument is less than 1
    """
    if M.map[0] < 1:
        raise ValueError("Argument of acosh must be >= 1")
        
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arccosh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i] / np.sqrt(M.map[0]**2 - 1)
    return result

def atanh(M):
    """
    Inverse hyperbolic tangent of a TPSAad object.
    
    Args:
        M: TPSAad object
    Returns:
        TPSAad object
    Raises:
        ValueError: if argument is outside (-1, 1)
    """
    if abs(M.map[0]) >= 1:
        raise ValueError("Argument of atanh must be in (-1, 1)")
        
    result = TPSAad()
    result.map = M.map.copy()
    result.map[0] = np.arctanh(M.map[0])
    for i in range(1, M.terms):
        result.map[i] = M.map[i] / (1 - M.map[0]**2)
    return result
# Add other mathematical functions (tan, sinh, cosh, etc.) similarly...
