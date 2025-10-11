#include <iostream>
#include <array>
#include "TPSAad.h" // Assuming this contains the TPSA AD implementation
		   
// Forward declaration of quadmap function
void quadmap(std::array<TPSAad, 4>& tpsaPtc, const TPSAad& tau, 
            const TPSAad& outin, const double xk);

int main() {
    // # of control variables
    const int noi = 7;

    // 7 control variables
    TPSAad x1, x2, x3, x4, x5, x6, x7;

    // 2D Twiss parameters (beta_x,alpha_x,beta_y,alpha_y)
    std::array<TPSAad, 4> twiss;

    // Initialize each lattice variable
    // drift length
    x1.assign(0.2, 1);
    // quad length
    x2.assign(0.1, 2);
    // quad strength (k)
    x3.assign(29.6, 3);
    // drift length
    x4.assign(0.4, 4);
    // quad length
    x5.assign(0.1, 5);
    // quad strength (k)
    x6.assign(-29.6, 6);
    // drift length
    x7.assign(0.2, 7);

    // initial Twiss parameters
    for(int i = 0; i < 4; ++i) {
        twiss[0] = 1.0;
        twiss[1] = 0.5;
        twiss[2] = 1.0;
        twiss[3] = -0.5;
    }

    // loop through nperd of the FODO lattice
    const int nperd = 1;
    for(int i = 0; i < nperd; ++i) {
        // drift
        double xk = 0.0;
        quadmap(twiss, x1, x2, xk);
        // quad
        xk = x3.map[0];
        quadmap(twiss, x2, x3, xk);
        // drift
        xk = 0.0;
        quadmap(twiss, x4, x2, xk);
        // quad
        xk = x6.map[0];
        quadmap(twiss, x5, x6, xk);
        // drift
        xk = 0.0;
        quadmap(twiss, x7, x2, xk);
    }

    std::cout << "Beta_x and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) std::cout << twiss[0].map[i] << " ";
    std::cout << std::endl;

    std::cout << "Alpha_x and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) std::cout << twiss[1].map[i] << " ";
    std::cout << std::endl;

    std::cout << "Beta_y and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) std::cout << twiss[2].map[i] << " ";
    std::cout << std::endl;

    std::cout << "Alpha_y and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) std::cout << twiss[3].map[i] << " ";
    std::cout << std::endl;

    return 0;
}

void quadmap(std::array<TPSAad, 4>& tpsaPtc, const TPSAad& tau, 
            const TPSAad& outin, const double xk) {
    std::array<TPSAad, 4> tpsaPtcTemp;
    TPSAad co, si, ch, sh, sqk, gammax, gammay;

    gammax = (1.0 + pow(tpsaPtc[1], 2)) / tpsaPtc[0];
    gammay = (1.0 + pow(tpsaPtc[3], 2)) / tpsaPtc[2];

    if (xk == 0.0) {
        tpsaPtc[0] = tpsaPtc[0] - 2*tpsaPtc[1]*tau + tau*tau*gammax;
        tpsaPtc[2] = tpsaPtc[2] - 2*tpsaPtc[3]*tau + tau*tau*gammay;
    }
    else if (xk > 0.0) {
        sqk = sqrt(outin);
        co = cos(sqk * tau);
        si = sin(sqk * tau);
        ch = cosh(sqk * tau);
        sh = sinh(sqk * tau);

        tpsaPtcTemp[0] = tpsaPtc[0]*co*co - 
                         2*tpsaPtc[1]*co*si/sqk + pow(si/sqk,2)*gammax;
        tpsaPtcTemp[1] = tpsaPtc[0]*co*si*sqk + 
                         tpsaPtc[1]*(co*co-si*si) - co*si/sqk*gammax;
        tpsaPtcTemp[2] = tpsaPtc[2]*ch*ch - 
                         2*tpsaPtc[3]*ch*sh/sqk + pow(sh/sqk,2)*gammay;
        tpsaPtcTemp[3] = -tpsaPtc[2]*ch*sh*sqk + 
                         tpsaPtc[3]*(ch*ch+sh*sh) - ch*sh/sqk*gammay;

        for(int i = 0; i < 4; ++i) {
            tpsaPtc[i] = tpsaPtcTemp[i];
        }
    }
    else if (xk < 0.0) {
        sqk = sqrt(-outin);
        co = cos(sqk * tau);
        si = sin(sqk * tau);
        ch = cosh(sqk * tau);
        sh = sinh(sqk * tau);

        tpsaPtcTemp[0] = tpsaPtc[0]*ch*ch - 
                         2*tpsaPtc[1]*ch*sh/sqk + pow(sh/sqk,2)*gammax;
        tpsaPtcTemp[1] = -tpsaPtc[0]*ch*sh*sqk + 
                         tpsaPtc[1]*(ch*ch+sh*sh) - ch*sh/sqk*gammax;
        tpsaPtcTemp[2] = tpsaPtc[2]*co*co - 
                         2*tpsaPtc[3]*co*si/sqk + pow(si/sqk,2)*gammay;
        tpsaPtcTemp[3] = tpsaPtc[2]*co*si*sqk + 
                         tpsaPtc[3]*(co*co-si*si) - co*si/sqk*gammay;

        for(int i = 0; i < 4; ++i) {
            tpsaPtc[i] = tpsaPtcTemp[i];
        }
    }
}
