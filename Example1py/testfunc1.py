# Created by Ji Qiang

from TPSAad import TPSAad, exp, sin, cos, sinh # Assuming class is implemented in a separate module
#import numpy as np

def main():
    # Create objects
#    x1 = TPSAad()
#    x2 = TPSAad()
#    x3 = TPSAad()
#    x4 = TPSAad()
#    x5 = TPSAad()
#    func = TPSAad()

    # Assign values and indices
#    x1.assign(3.0, 1)
#    x2.assign(0.1, 2)
#    x3.assign(0.1, 3)
#    x4.assign(1.0, 4)
#    x5.assign(1.0, 5)
     
    x1 = TPSAad(3.0, 1)
    x2 = TPSAad(0.1, 2)
    x3 = TPSAad(0.1, 3)
    x4 = TPSAad(1.0, 4)
    x5 = TPSAad(1.0, 5)

    # Test function
    func = (2*cos(x1/x2) + x1/x2 + exp(x2))*x3 + 2*x4 + sinh(x5)

    # Print results
    print("func. value and its derivatives w.r.t. 5 variables")
    print(f"{func.map[0]} {func.map[1]} {func.map[2]} {func.map[3]} {func.map[4]} {func.map[5]}")

if __name__ == "__main__":
    main()
