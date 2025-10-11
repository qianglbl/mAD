//****************************
//
//*** Copyright Notice ***
//
//A multi-language auto differentiation package (mAD) Copyright (c) 2025, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
//
//If you have questions about your rights to use or distribute this software,
//please contact Berkeley Lab's Intellectual Property Office at
//IPO@lbl.gov.
//
//NOTICE.  This Software was developed under funding from the U.S. Department
//of Energy and the U.S. Government consequently retains certain rights.  As
//such, the U.S. Government has been granted for itself and others acting on
//its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
//Software to reproduce, distribute copies to the public, prepare derivative
//works, and perform publicly and display publicly, and to permit others to do so.
//
//****************************
// Java class for forward auto-differentiation
// Ji Qiang with help of Antropic Claude Sonnet - LBNL, jqiang@lbl.gov

public class TPSAad {
    //dimmax: maximum # of variables to be differentiable
    //For better performance, this number should be the # differentiable variables.
    public static final int DIMMAX = 7;
    public double[] map;
    public int terms;
    //private int terms;

    // Constructors
    public TPSAad() {
        terms = DIMMAX + 1;
        map = new double[terms];
    }

    public TPSAad(int nvar) {
        terms = nvar + 1;
        map = new double[terms];
    }

    public TPSAad(double a) {
        terms = DIMMAX + 1;
        map = new double[terms];
        map[0] = a;
    }

    public TPSAad(double a, int ivar) {
        terms = DIMMAX + 1;
        map = new double[terms];
        map[0] = a;
        map[ivar] = 1.0;
    }

    public TPSAad copy() {
        TPSAad result = new TPSAad();
        System.arraycopy(this.map, 0, result.map, 0, this.terms);
        return result;
    }

    // Copy constructor
    public TPSAad(TPSAad other) {
        this.terms = other.terms;
        this.map = new double[terms];
        System.arraycopy(other.map, 0, this.map, 0, terms);
    }

    // Getters
    public int getTerms() {
        return terms;
    }

    public int getNvar() {
        return terms - 1;
    }

    // Utility methods
    public double getValue() {
        return map[0];
    }

    // Assignment methods
    public void assign(double a) {
        terms = DIMMAX + 1;
        map[0] = a;
        for (int i = 1; i < terms; i++) {
            map[i] = 0.0;
        }
    }

    public void assign(double a, int ivar) {
        terms = DIMMAX + 1;
        for (int i = 0; i < terms; i++) {
            map[i] = 0.0;
        }
        map[ivar] = 1.0;
        map[0] = a;
    }

    // Basic arithmetic operations
    public TPSAad add(TPSAad other) {
        TPSAad result = new TPSAad(this.terms - 1);
        for (int i = 0; i < terms; i++) {
            result.map[i] = this.map[i] + other.map[i];
        }
        return result;
    }

    public static TPSAad add(double a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Add scalar to constant term
        result.map[0] = a + b.map[0];
        // Copy derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = b.map[i];
        }
        return result;
    }

    public static TPSAad add(int a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Add scalar to constant term
        result.map[0] = a + b.map[0];
        // Copy derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = b.map[i];
        }
        return result;
    }

    public static TPSAad add(TPSAad a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Add corresponding terms including derivatives
        for (int i = 0; i < b.terms; i++) {
            result.map[i] = a.map[i] + b.map[i];
        }
        return result;
    }

    public static TPSAad add(TPSAad a, double b) {
        return add(b, a); // Reuse the double + TPSAad implementation
    }

    public static TPSAad add(TPSAad a, int b) {
        return add(b, a); // Reuse the double + TPSAad implementation
    }

    public void addInPlace(double value) {
        this.map[0] += value;
    }

    public void addInPlace(TPSAad other) {
        for (int i = 0; i < this.terms; i++) {
            this.map[i] += other.map[i];
        }
    }


    public TPSAad subtract(TPSAad other) {
        TPSAad result = new TPSAad(this.terms - 1);
        for (int i = 0; i < terms; i++) {
            result.map[i] = this.map[i] - other.map[i];
        }
        return result;
    }

    public static TPSAad subtract(TPSAad a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Subtract corresponding terms including derivatives
        for (int i = 0; i < a.terms; i++) {
            result.map[i] = a.map[i] - b.map[i];
        }
        return result;
    }

    public static TPSAad subtract(TPSAad a, double b) {
        TPSAad result = new TPSAad();
        // Subtract scalar from constant term
        result.map[0] = a.map[0] - b;
        // Copy derivatives from a (derivatives of a constant are 0)
        for (int i = 1; i < a.terms; i++) {
            result.map[i] = a.map[i];
        }
        return result;
    }

    public static TPSAad subtract(TPSAad a, int b) {
        TPSAad result = new TPSAad();
        // Subtract scalar from constant term
        result.map[0] = a.map[0] - b;
        // Copy derivatives from a (derivatives of a constant are 0)
        for (int i = 1; i < a.terms; i++) {
            result.map[i] = a.map[i];
        }
        return result;
    }

    public static TPSAad subtract(double a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Subtract constant term
        result.map[0] = a - b.map[0];
        // Negate derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = -b.map[i];
        }
        return result;
    }

    public static TPSAad subtract(int a, TPSAad b) {
        TPSAad result = new TPSAad();
        // Subtract constant term
        result.map[0] = a - b.map[0];
        // Negate derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = -b.map[i];
        }
        return result;
    }

/*
    public static TPSAad subtract(TPSAad b, double a) {
        TPSAad result = new TPSAad();
        // Subtract constant term
        result.map[0] = b.map[0] - a;
        // Negate derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = b.map[i];
        }
        return result;
    }

    public static TPSAad subtract(TPSAad b, int a) {
        TPSAad result = new TPSAad();
        // Subtract constant term
        result.map[0] = b.map[0] - a;
        // Negate derivatives from b
        for (int i = 1; i < b.terms; i++) {
            result.map[i] = b.map[i];
        }
        return result;
    }
*/

    public void subtractInPlace(TPSAad other) {
        for (int i = 0; i < this.terms; i++) {
            this.map[i] -= other.map[i];
        }
    }

    public void subtractInPlace(double value) {
        this.map[0] -= value;
    }
 
    public TPSAad negate() {
        TPSAad result = new TPSAad();
        for (int i = 0; i < this.terms; i++) {
            result.map[i] = -this.map[i];
        }
        return result;
    }

    public TPSAad multiply(TPSAad other) {
        TPSAad result = new TPSAad(this.terms - 1);
        result.map[0] = this.map[0] * other.map[0];
        for (int i = 1; i < terms; i++) {
            result.map[i] = other.map[0] * this.map[i] + this.map[0] * other.map[i];
        }
        return result;
    }

   public static TPSAad multiply(TPSAad a, TPSAad b) {
        TPSAad result = new TPSAad();

        // Constant term (zeroth order)
        result.map[0] = a.map[0] * b.map[0];

        // First order derivatives using product rule
        // d(u*v) = u*dv + v*du
        for (int i = 1; i < a.terms; i++) {
            result.map[i] = a.map[0] * b.map[i] + b.map[0] * a.map[i];
        }

        return result;
    }

    public static TPSAad multiply(TPSAad a, double b) {
        TPSAad result = new TPSAad();
        // Multiply all terms by scalar
        for (int i = 0; i < a.terms; i++) {
            result.map[i] = a.map[i] * b;
        }
        return result;
    }

    public static TPSAad multiply(double a, TPSAad b) {
        return multiply(b, a); // Reuse TPSAad * scalar implementation
    }

    public static TPSAad multiply(TPSAad a, int b) {
        TPSAad result = new TPSAad();
        // Multiply all terms by scalar
        for (int i = 0; i < a.terms; i++) {
            result.map[i] = a.map[i] * b;
        }
        return result;
    }

    public static TPSAad multiply(int a, TPSAad b) {
        return multiply(b, a); // Reuse TPSAad * scalar implementation
    }


    public void multiplyInPlace(TPSAad other) {
        double[] tempMap = new double[this.terms];

        // Calculate all products
        tempMap[0] = this.map[0] * other.map[0];

        for (int i = 1; i < this.terms; i++) {
            tempMap[i] = this.map[0] * other.map[i] +
                        other.map[0] * this.map[i];
        }

        // Copy results back to this object
        System.arraycopy(tempMap, 0, this.map, 0, this.terms);
    }

    public void multiplyInPlace(double scalar) {
        for (int i = 0; i < this.terms; i++) {
            this.map[i] *= scalar;
        }
    }

    public void multiplyInPlace(int scalar) {
        for (int i = 0; i < this.terms; i++) {
            this.map[i] *= scalar;
        }
    }

    public TPSAad square() {
        return multiply(this, this);
    }

    public TPSAad pow(int n) {
        if (n == 0) {
            TPSAad result = new TPSAad();
            result.map[0] = 1.0;
            return result;
        }
        if (n == 1) {
            return this.copy();
        }
        if (n == 2) {
            return this.square();
        }

        // For higher powers, use recursive multiplication
        TPSAad result = this.copy();
        for (int i = 1; i < n; i++) {
            result = multiply(result, this);
        }
        return result;
    }

    public TPSAad divide(TPSAad other) {
        if (Math.abs(other.map[0]) < 1.0e-14) {
            throw new ArithmeticException("Division by zero in TPSAad");
        }
        TPSAad result = new TPSAad(this.terms - 1);
        result.map[0] = this.map[0] / other.map[0];
        for (int i = 1; i < terms; i++) {
            result.map[i] = (this.map[i] * other.map[0] - other.map[i] * this.map[0]) 
                           / (other.map[0] * other.map[0]);
        }
        return result;
    }

    public static TPSAad divide(TPSAad a, TPSAad b) {
        if (Math.abs(b.map[0]) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        TPSAad result = new TPSAad();

        // Constant term
        result.map[0] = a.map[0] / b.map[0];

        // First order derivatives using quotient rule
        // d(u/v) = (v*du - u*dv)/(v*v)
        double v = b.map[0];
        double v2 = v * v;

        for (int i = 1; i < a.terms; i++) {
            result.map[i] = (v * a.map[i] - a.map[0] * b.map[i]) / v2;
        }

        return result;
    }

    /**
     * Divides a TPSAad object by a scalar
     * @param a TPSAad object (numerator)
     * @param b scalar value (denominator)
     * @return new TPSAad object containing the quotient
     * @throws ArithmeticException if dividing by zero
     */
    public static TPSAad divide(TPSAad a, double b) {
        if (Math.abs(b) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        TPSAad result = new TPSAad();
        // Divide all terms by scalar
        for (int i = 0; i < a.terms; i++) {
            result.map[i] = a.map[i] / b;
        }
        return result;
    }

    public static TPSAad divide(TPSAad a, int b) {
        if (Math.abs(b) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        TPSAad result = new TPSAad();
        // Divide all terms by scalar
        for (int i = 0; i < a.terms; i++) {
            result.map[i] = a.map[i] / b;
        }
        return result;
    }

    /**
     * Divides a scalar by a TPSAad object
     * @param a scalar value (numerator)
     * @param b TPSAad object (denominator)
     * @return new TPSAad object containing the quotient
     * @throws ArithmeticException if dividing by zero
     */
    public static TPSAad divide(double a, TPSAad b) {
        if (Math.abs(b.map[0]) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        TPSAad result = new TPSAad();

        // Constant term
        result.map[0] = a / b.map[0];

        // Derivatives using chain rule
        // d(a/v) = -a*dv/v^2
        double v = b.map[0];
        double v2 = v * v;

        for (int i = 1; i < b.terms; i++) {
            result.map[i] = -a * b.map[i] / v2;
        }

        return result;
    }

    public static TPSAad divide(int a, TPSAad b) {
        if (Math.abs(b.map[0]) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        TPSAad result = new TPSAad();

        // Constant term
        result.map[0] = a / b.map[0];

        // Derivatives using chain rule
        // d(a/v) = -a*dv/v^2
        double v = b.map[0];
        double v2 = v * v;

        for (int i = 1; i < b.terms; i++) {
            result.map[i] = -a * b.map[i] / v2;
        }

        return result;
    }

    /**
     * Divides this TPSAad object by another in place
     * @param other TPSAad object to divide by
     * @throws ArithmeticException if dividing by zero
     */
    public void divideInPlace(TPSAad other) {
        if (Math.abs(other.map[0]) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        double[] tempMap = new double[this.terms];

        // Calculate quotient
        double v = other.map[0];
        double v2 = v * v;

        tempMap[0] = this.map[0] / v;

        for (int i = 1; i < this.terms; i++) {
            tempMap[i] = (v * this.map[i] - this.map[0] * other.map[i]) / v2;
        }

        // Copy results back to this object
        System.arraycopy(tempMap, 0, this.map, 0, this.terms);
    }

    /**
     * Divides this TPSAad object by a scalar in place
     * @param scalar value to divide by
     * @throws ArithmeticException if dividing by zero
     */
    public void divideInPlace(double scalar) {
        if (Math.abs(scalar) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        for (int i = 0; i < this.terms; i++) {
            this.map[i] /= scalar;
        }
    }

    public void divideInPlace(int scalar) {
        if (Math.abs(scalar) < 1.0e-14) {
            throw new ArithmeticException("Division by zero or near-zero value");
        }

        for (int i = 0; i < this.terms; i++) {
            this.map[i] /= scalar;
        }
    }

    /**
     * Computes the reciprocal (1/x) of this TPSAad object
     * @return new TPSAad object containing the reciprocal
     * @throws ArithmeticException if this value is zero
     */
    public TPSAad reciprocal() {
        return divide(1.0, this);
    }

    // Static math functions
    public static TPSAad exp(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        double temp = M.map[0];
        result.map[0] = Math.exp(temp);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = Math.exp(M.map[0]) * M.map[i];
        }
        return result;
    }

    public static TPSAad sin(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.sin(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = Math.cos(M.map[0]) * M.map[i];
        }
        return result;
    }

    public static TPSAad cos(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.cos(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = -Math.sin(M.map[0]) * M.map[i];
        }
        return result;
    }

    // Additional mathematical functions
    public static TPSAad tan(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.tan(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = (1.0 + Math.pow(Math.tan(M.map[0]), 2)) * M.map[i];
        }
        return result;
    }

    public static TPSAad sinh(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.sinh(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = Math.cosh(M.map[0]) * M.map[i];
        }
        return result;
    }

    public static TPSAad cosh(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.cosh(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = Math.sinh(M.map[0]) * M.map[i];
        }
        return result;
    }

    public static TPSAad tanh(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.tanh(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = (1.0 - Math.pow(Math.tanh(M.map[0]), 2)) * M.map[i];
        }
        return result;
    }

    public static TPSAad asin(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.asin(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = M.map[i] / Math.sqrt(1 - Math.pow(M.map[0], 2));
        }
        return result;
    }
    public static TPSAad acos(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.acos(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = -M.map[i] / Math.sqrt(1 - Math.pow(M.map[0], 2));
        }
        return result;
    }

    public static TPSAad atan(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.atan(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = M.map[i] / (1 + Math.pow(M.map[0], 2));
        }
        return result;
    }

    public static TPSAad pow(TPSAad M, double a) {
        TPSAad result = new TPSAad(M.terms - 1);
        if (Math.abs(a - 1.0) < 1e-14) {
            return new TPSAad(M);
        } else if (Math.abs(a) < 1e-14) {
            result.map[0] = 1.0;
            return result;
        }

        result.map[0] = Math.pow(M.map[0], a);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = a * Math.pow(M.map[0], a - 1) * M.map[i];
        }
        return result;
    }

    public static TPSAad pow(TPSAad M, int a) {
        TPSAad result = new TPSAad(M.terms - 1);
        if (a == 1) {
            return new TPSAad(M);
        } else if (a == 0) {
            result.map[0] = 1.0;
            return result;
        }

        result.map[0] = Math.pow(M.map[0], a);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = a * Math.pow(M.map[0], a - 1) * M.map[i];
        }
        return result;
    }

    public static TPSAad sqrt(TPSAad M) {
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.sqrt(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = 0.5 * M.map[i] / Math.sqrt(M.map[0]);
        }
        return result;
    }

    public static TPSAad log(TPSAad M) {
        if (Math.abs(M.map[0]) < 1.0e-15) {
            throw new ArithmeticException("Log of zero in TPSAad");
        }
        TPSAad result = new TPSAad(M.terms - 1);
        result.map[0] = Math.log(M.map[0]);
        for (int i = 1; i < M.terms; i++) {
            result.map[i] = M.map[i] / M.map[0];
        }
        return result;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Value: ").append(map[0]).append("\n");
        sb.append("Derivatives: [");
        for (int i = 1; i < terms; i++) {
            sb.append(map[i]);
            if (i < terms - 1) sb.append(", ");
        }
        sb.append("]");
        return sb.toString();
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof TPSAad)) return false;

        TPSAad other = (TPSAad) obj;
        if (this.terms != other.terms) return false;

        for (int i = 0; i < terms; i++) {
            if (Math.abs(this.map[i] - other.map[i]) > 1e-14) return false;
        }
        return true;
    }

}
