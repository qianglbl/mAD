from TPSAad import * # Assuming class is implemented i n a separate module
from typing import List

def main():
    # Number of control variables
    noi = 7
    
    # Initialize control variables
    #x1, x2, x3, x4, x5, x6, x7 = [TPSAad() for _ in range(7)]

    # Initialize control variables
    x1 = TPSAad(0.2, 1)  # drift length
    x2 = TPSAad(0.1, 2)  # quad length
    x3 = TPSAad(29.6, 3)  # quad strength (k)
    x4 = TPSAad(0.4, 4)  # drift length
    x5 = TPSAad(0.1, 5)  # quad length
    x6 = TPSAad(-29.6, 6)  # quad strength (k)
    x7 = TPSAad(0.2, 7)  # drift length

    
    # Initialize Twiss parameters
    twiss = [TPSAad() for _ in range(4)]
    
    # Initial Twiss parameters
    twiss[0]=TPSAad(1.0)
    twiss[1]=TPSAad(0.5)
    twiss[2]=TPSAad(1.0)
    twiss[3]=TPSAad(-0.5)
    
    # Loop through FODO lattice
    nperd = 1
    for _ in range(nperd):
        # drift
        xk = 0.0
        quadmap(twiss, x1, x2, xk)
        # quad
        xk = x3.map[0]
        quadmap(twiss, x2, x3, xk)
        # drift
        xk = 0.0
        quadmap(twiss, x4, x2, xk)
        # quad
        xk = x6.map[0]
        quadmap(twiss, x5, x6, xk)
        # drift
        xk = 0.0
        quadmap(twiss, x7, x2, xk)
    
    # Print results
    coord_names = ['BetaX', 'AlphaX', 'BetaY', 'AlphaY']
    for i, name in enumerate(coord_names):
        print(f"{name} coordinate and its derivatives w.r.t. 7 variables:")
        print(' '.join(f"{twiss[i].map[j]}" for j in range(noi + 1)))
        print()


def quadmap(tpsaPtc: List[TPSAad], tau: TPSAad, outin: TPSAad, xk: float) -> None:
    """Implement quad map transformation"""
    tpsaPtcTemp = [TPSAad() for _ in range(4)]
    
    # Calculate gamma
    gammax = (1.0 + pow(tpsaPtc[1],2))/tpsaPtc[0]
    gammay = (1.0 + pow(tpsaPtc[3],2))/tpsaPtc[2]
    
    if abs(xk) < 1e-10:
        tpsaPtc[0] = tpsaPtc[0] - 2*tpsaPtc[1]*tau + tau*tau*gammax
        tpsaPtc[2] = tpsaPtc[2] - 2*tpsaPtc[3]*tau + tau*tau*gammay
    
    elif xk > 0.0:
        sqk = sqrt(outin)
        co = cos(sqk * tau)
        si = sin(sqk * tau)
        ch = cosh(sqk * tau)
        sh = sinh(sqk * tau)
        
        # Update temp values
        tpsaPtcTemp[0] = tpsaPtc[0]*(co*co) - 2*tpsaPtc[1]*(co*si/sqk) + pow((si/sqk),2)*gammax
        tpsaPtcTemp[1] = tpsaPtc[0]*(co*si*sqk) + tpsaPtc[1]*(co*co-si*si) - (co*si/sqk)*gammax
        tpsaPtcTemp[2] = tpsaPtc[2]*(ch*ch) - 2*tpsaPtc[3]*(ch*sh/sqk) + pow((sh/sqk),2)*gammay
        tpsaPtcTemp[3] = -tpsaPtc[2]*(ch*sh*sqk) + tpsaPtc[3]*(ch*ch+sh*sh) - (ch*sh/sqk)*gammay
        
        # Update original array
        for i in range(4):
            tpsaPtc[i] = tpsaPtcTemp[i]
    
    elif xk < 0.0:
        sqk = sqrt(-outin)
        co = cos(sqk * tau)
        si = sin(sqk * tau)
        ch = cosh(sqk * tau)
        sh = sinh(sqk * tau)
        
        # Update temp values
        tpsaPtcTemp[0] = tpsaPtc[0]*(ch*ch) - 2*tpsaPtc[1]*(ch*sh/sqk) + pow((sh/sqk),2)*gammax
        tpsaPtcTemp[1] = -tpsaPtc[0]*(ch*sh*sqk) + tpsaPtc[1]*(ch*ch+sh*sh) - (ch*sh/sqk)*gammax
        tpsaPtcTemp[2] = tpsaPtc[2]*(co*co) - 2*tpsaPtc[3]*(co*si/sqk) + pow((si/sqk),2)*gammay
        tpsaPtcTemp[3] = tpsaPtc[2]*(co*si*sqk) + tpsaPtc[3]*(co*co-si*si) - (co*si/sqk)*gammay
        
        # Update original array
        for i in range(4):
            tpsaPtc[i] = tpsaPtcTemp[i]


if __name__ == "__main__":
    main()
