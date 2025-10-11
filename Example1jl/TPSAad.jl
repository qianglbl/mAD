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
# Julia module for forward auto-differentiation
# Ji Qiang with help of Antropic Claude Sonnet - LBNL, jqiang@lbl.gov
#

module TPSAad

export DTPSAD
export dtpsad_initialize,init_tpsad,assign_v1!,assign_v2!
export +, -, *, /, ==, ^
export exp, log, sqrt, sin, cos, tan
export asin, acos, atan, sinh, cosh, tanh
export asinh, acosh, atanh

#dimmax: maximum # of variables to be differentiable
#For better performance, this number should be the # differentiable variables.
const dimmax = 5
global TPS_Dim = 0

mutable struct DTPSAD
    terms::Int32
    map::Vector{Float64}
    
    function DTPSAD()
        new(dimmax+1, zeros(Float64, dimmax + 1))
    end
end

# Helper functions for initialization
function dtpsad_initialize(nvar::Integer)
    global TPS_Dim = nvar
end

function init_tpsad(this::DTPSAD, nvar::Integer)
    this.terms = nvar + 1
    this.map[1:this.terms] .= 0.0
end

# Assignment operations
function assign_v1!(this::DTPSAD, a::Float64)
    this.terms = TPS_Dim + 1
    this.map[1] = a
    this.map[2:this.terms] .= 0.0
end

function assign_v2!(this::DTPSAD, a::Float64, i_var::Integer)
    this.terms = TPS_Dim + 1
    this.map[1] = a
    this.map[2:this.terms] .= 0.0
    this.map[i_var + 1] = 1.0
end

# Operator overloading
import Base: +, -, *, /, ==, ^
import Base: exp, log, sqrt, sin, cos, tan
import Base: asin, acos, atan, sinh, cosh, tanh
import Base: asinh, acosh, atanh

# Addition operators
function +(M::DTPSAD, N::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = M.map[1:M.terms] + N.map[1:M.terms]
    return result
end

function +(M::DTPSAD, N::Number)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = M.map[1] + N
    result.map[2:M.terms] = M.map[2:M.terms]
    return result
end

+(N::Number, M::DTPSAD) = M + N

function +(M::DTPSAD)  # Unary plus
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = M.map[1:M.terms]
    return result
end

# Subtraction operators
function -(M::DTPSAD, N::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = M.map[1:M.terms] - N.map[1:M.terms]
    return result
end

function -(M::DTPSAD, N::Number)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = M.map[1] - N
    result.map[2:M.terms] = M.map[2:M.terms]
    return result
end

function -(N::Number, M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = N - M.map[1]
    result.map[2:M.terms] = -M.map[2:M.terms]
    return result
end

function -(M::DTPSAD)  # Unary minus
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = -M.map[1:M.terms]
    return result
end

# Multiplication operators
function *(M::DTPSAD, N::DTPSAD)
    result = DTPSAD()
    result.terms = max(M.terms, N.terms)
    result.map[1] = M.map[1] * N.map[1]
    result.map[2:M.terms] = M.map[1] * N.map[2:N.terms] +
                           M.map[2:M.terms] * N.map[1]
    return result
end

function *(M::DTPSAD, N::Number)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = M.map[1:M.terms] * N
    return result
end

*(N::Number, M::DTPSAD) = M * N

# Division operators
function /(M::DTPSAD, N::DTPSAD)
    if abs(N.map[1]) < 1.0e-14
        error("Error: Divide by zero in DTPSAD division")
    end

    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = M.map[1] / N.map[1]
    result.map[2:M.terms] = (M.map[2:M.terms] * N.map[1] -
                            N.map[2:M.terms] * M.map[1]) / (N.map[1]^2)
    return result
end

function /(M::DTPSAD, N::Number)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1:M.terms] = M.map[1:M.terms] / N
    return result
end

function /(N::Number, M::DTPSAD)
    if abs(M.map[1]) < 1.0e-14
        error("Error: Divide by zero in DTPSAD division")
    end

    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = N / M.map[1]
    result.map[2:M.terms] = -M.map[2:M.terms] * N / (M.map[1]^2)
    return result
end

# Comparison
function ==(M::DTPSAD, N::DTPSAD)
    if M.terms != N.terms
        return false
    end

    for i in 1:M.terms
        if abs(M.map[i]) > 1.0e-30 || abs(N.map[i]) > 1.0e-30
            if abs(M.map[i] / N.map[i] - 1.0) > 1.0e-15
                return false
            end
        end
    end
    return true
end

# Inverse function
function inv(M::DTPSAD)
    if abs(M.map[1]) < 1.0e-14
        error("Error: Divide by zero in DTPSAD inverse")
    end

    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = 1.0 / M.map[1]
    result.map[2:M.terms] = -M.map[2:M.terms] / (M.map[1]^2)
    return result
end

# Power functions
function ^(M::DTPSAD, b::Float64)
    result = DTPSAD()
    result.terms = M.terms

    if abs(b - 1.0) < 1.0e-14
        return M
    elseif abs(b) < 1.0e-14
        result.map[1] = 1.0
        result.map[2:M.terms] .= 0.0
        return result
    else
        result.map[1] = M.map[1]^b
        result.map[2:M.terms] = b * (M.map[1]^(b-1)) * M.map[2:M.terms]
        return result
    end
end

function ^(M::DTPSAD, b::Integer)
    result = DTPSAD()
    result.terms = M.terms

    if b == 1
        return M
    elseif b == 0
        result.map[1] = 1.0
        result.map[2:M.terms] .= 0.0
        return result
    else
        result.map[1] = M.map[1]^b
        result.map[2:M.terms] = b * (M.map[1]^(b-1)) * M.map[2:M.terms]
        return result
    end
end

# Exponential function
function exp(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = exp(M.map[1])
    result.map[2:M.terms] = exp(M.map[1]) * M.map[2:M.terms]
    return result
end

# Natural logarithm
function log(M::DTPSAD)
    if abs(M.map[1]) < 1.0e-15
        error("Error: Log of zero in DTPSAD")
    end
    
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = log(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / M.map[1]
    return result
end

# Square root
function sqrt(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = sqrt(M.map[1])
    result.map[2:M.terms] = 0.5 * M.map[2:M.terms] / sqrt(M.map[1])
    return result
end

# Trigonometric functions
function sin(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    sin_a0 = sin(M.map[1])
    cos_a0 = cos(M.map[1])
    result.map[1] = sin_a0
    result.map[2:M.terms] = cos_a0 * M.map[2:M.terms]
    return result
end

function cos(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    sin_a0 = sin(M.map[1])
    cos_a0 = cos(M.map[1])
    result.map[1] = cos_a0
    result.map[2:M.terms] = -sin_a0 * M.map[2:M.terms]
    return result
end

function tan(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = tan(M.map[1])
    result.map[2:M.terms] = (1.0 + tan(M.map[1])^2) * M.map[2:M.terms]
    return result
end

# Inverse trigonometric functions
function asin(M::DTPSAD)
    if abs(M.map[1]) >= 1
        error("Error: Invalid argument for asin in DTPSAD")
    end
    
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = asin(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / sqrt(1 - M.map[1]^2)
    return result
end

function acos(M::DTPSAD)
    if abs(M.map[1]) >= 1
        error("Error: Invalid argument for acos in DTPSAD")
    end
    
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = acos(M.map[1])
    result.map[2:M.terms] = -M.map[2:M.terms] / sqrt(1 - M.map[1]^2)
    return result
end

function atan(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = atan(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / (1 + M.map[1]^2)
    return result
end

# Hyperbolic functions
function sinh(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = sinh(M.map[1])
    result.map[2:M.terms] = cosh(M.map[1]) * M.map[2:M.terms]
    return result
end

function cosh(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = cosh(M.map[1])
    result.map[2:M.terms] = sinh(M.map[1]) * M.map[2:M.terms]
    return result
end

function tanh(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = tanh(M.map[1])
    result.map[2:M.terms] = (1.0 - tanh(M.map[1])^2) * M.map[2:M.terms]
    return result
end

# Inverse hyperbolic functions
function asinh(M::DTPSAD)
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = asinh(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / sqrt(1 + M.map[1]^2)
    return result
end

function acosh(M::DTPSAD)
    if M.map[1] <= 1
        error("Error: Invalid argument for acosh in DTPSAD")
    end
    
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = acosh(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / sqrt(M.map[1]^2 - 1)
    return result
end

function atanh(M::DTPSAD)
    if abs(M.map[1]) >= 1
        error("Error: Invalid argument for atanh in DTPSAD")
    end
    
    result = DTPSAD()
    result.terms = M.terms
    result.map[1] = atanh(M.map[1])
    result.map[2:M.terms] = M.map[2:M.terms] / (1 - M.map[1]^2)
    return result
end

end 
