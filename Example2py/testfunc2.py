from TPSAad import * # Assuming class is implemented i n a separate module
from typing import List

def main():
    # # of control variables
    noi = 7
    
    # Initialize control variables
    x1 = TPSAad(0.2, 1)  # drift length
    x2 = TPSAad(0.1, 2)  # quad length
    x3 = TPSAad(29.6, 3)  # quad strength (k)
    x4 = TPSAad(0.4, 4)  # drift length
    x5 = TPSAad(0.1, 5)  # quad length
    x6 = TPSAad(-29.6, 6)  # quad strength (k)
    x7 = TPSAad(0.2, 7)  # drift length
    
    # Initialize particle coordinates (x,px,y,py)
    pt = [TPSAad() for _ in range(4)]
    pt[0] = TPSAad(1.0e-3)
    pt[1] = TPSAad(1.0e-3)
    pt[2] = TPSAad(-1.0e-3)
    pt[3] = TPSAad(-1.0e-3)
    
    # Loop through nperd of the FODO lattice
    nperd = 1
    for _ in range(nperd):
        # drift
        xk = 0.0
        quadmap(pt, x1, x2, xk)
        # quad
        xk = x3.map[0]
        quadmap(pt, x2, x3, xk)
        # drift
        xk = 0.0
        quadmap(pt, x4, x2, xk)
        # quad
        xk = x6.map[0]
        quadmap(pt, x5, x6, xk)
        # drift
        xk = 0.0
        quadmap(pt, x7, x2, xk)
    
    # Print results
    coord_names = ['X', 'Px', 'Y', 'Py']
    for i, name in enumerate(coord_names):
        print(f"{name} coordinate and its derivatives w.r.t. 7 variables:")
        print(' '.join(f"{pt[i].map[j]}" for j in range(noi + 1)))
        print()

def quadmap(tpsaPtc: List[TPSAad], tau: TPSAad, outin: TPSAad, xk: float) -> None:
    """
    Quadrupole map function for beam dynamics calculations.
    
    Args:
        tpsaPtc: List of 4 TPSAad objects representing particle coordinates
        tau: TPSAad object representing time step
        outin: TPSAad object
        xk: float representing quadrupole strength
    """
    tpsaPtcTemp = [TPSAad() for _ in range(4)]
    
    if abs(xk) < 1e-10:  # xk == 0.0
        tpsaPtc[0] = tpsaPtc[0] + tpsaPtc[1] * tau
        tpsaPtc[2] = tpsaPtc[2] + tpsaPtc[3] * tau
    
    elif xk > 0.0:
        sqk = sqrt(outin)
        co = cos(sqk * tau)
        si = sin(sqk * tau)
        ch = cosh(sqk * tau)
        sh = sinh(sqk * tau)
        
        tpsaPtcTemp[0] = tpsaPtc[0] * co + tpsaPtc[1] * si / sqk
        tpsaPtcTemp[1] = -tpsaPtc[0] * si * sqk + tpsaPtc[1] * co
        tpsaPtcTemp[2] = tpsaPtc[2] * ch + tpsaPtc[3] * sh / sqk
        tpsaPtcTemp[3] = tpsaPtc[2] * sh * sqk + tpsaPtc[3] * ch
        
        for i in range(4):
            tpsaPtc[i] = tpsaPtcTemp[i]
    
    elif xk < 0.0:
        sqk = sqrt(-outin)
        co = cos(sqk * tau)
        si = sin(sqk * tau)
        ch = cosh(sqk * tau)
        sh = sinh(sqk * tau)
        
        tpsaPtcTemp[0] = tpsaPtc[0] * ch + tpsaPtc[1] * sh / sqk
        tpsaPtcTemp[1] = tpsaPtc[0] * sh * sqk + tpsaPtc[1] * ch
        tpsaPtcTemp[2] = tpsaPtc[2] * co + tpsaPtc[3] * si / sqk
        tpsaPtcTemp[3] = -tpsaPtc[2] * si * sqk + tpsaPtc[3] * co
        
        for i in range(4):
            tpsaPtc[i] = tpsaPtcTemp[i]


if __name__ == "__main__":
    main()

