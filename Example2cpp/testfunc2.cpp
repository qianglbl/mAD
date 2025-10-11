//
// Created by Ji Qiang.
//
#include "TPSAad.h"
#include <iostream>
#include <array>

void quadmap(std::array<TPSAad, 4>& tpsaPtc, const TPSAad& tau,
            const TPSAad& outin, const double xk);

int main() {
    // # of control variables
    const int noi = 7;

    // 7 control variables
    TPSAad x1, x2, x3, x4, x5, x6, x7;

    // particle coordinates (x,px,y,py)
    std::array<TPSAad, 4> pt;

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

    // initial particle coordinates
    for(int i = 0; i < 4; ++i) {
        pt[0] = 1.0e-3;
        pt[1] = 1.0e-3;
        pt[2] = -1.0e-3;
        pt[3] = -1.0e-3;
    }

    // loop through nperd of the FODO lattice
    const int nperd = 1;
    for(int i = 0; i < nperd; ++i) {
        // drift
        double xk = 0.0;
        quadmap(pt, x1, x2, xk);
        // quad
        xk = x3.map[0];
        quadmap(pt, x2, x3, xk);
        // drift
        xk = 0.0;
        quadmap(pt, x4, x2, xk);
        // quad
        xk = x6.map[0];
        quadmap(pt, x5, x6, xk);
        // drift
        xk = 0.0;
        quadmap(pt, x7, x2, xk);
    }

    std::cout << "X coordinate and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) {
        std::cout << pt[0].map[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Px coordinate and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) {
        std::cout << pt[1].map[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Y coordinate and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) {
        std::cout << pt[2].map[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Py coordinate and its derivatives w.r.t. 7 variables:" << std::endl;
    for(int i = 0; i < noi+1; ++i) {
        std::cout << pt[3].map[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}

void quadmap(std::array<TPSAad, 4>& tpsaPtc, const TPSAad& tau, 
            const TPSAad& outin, const double xk) {
    
    std::array<TPSAad, 4> tpsaPtcTemp;
    TPSAad co, si, ch, sh, sqk;

    if (xk == 0.0) {
        tpsaPtc[0] = tpsaPtc[0] + tpsaPtc[1] * tau;
        tpsaPtc[2] = tpsaPtc[2] + tpsaPtc[3] * tau;
    }
    else if (xk > 0.0) {
        sqk = sqrt(outin);
        co = cos(sqk * tau);
        si = sin(sqk * tau);
        ch = cosh(sqk * tau);
        sh = sinh(sqk * tau);

        tpsaPtcTemp[0] = tpsaPtc[0] * co + 
                         tpsaPtc[1] * si / sqk;
        tpsaPtcTemp[1] = -tpsaPtc[0] * si * sqk + 
                         tpsaPtc[1] * co;
        tpsaPtcTemp[2] = tpsaPtc[2] * ch + 
                         tpsaPtc[3] * sh / sqk;
        tpsaPtcTemp[3] = tpsaPtc[2] * sh * sqk + 
                         tpsaPtc[3] * ch;

        tpsaPtc[0] = tpsaPtcTemp[0];
        tpsaPtc[1] = tpsaPtcTemp[1];
        tpsaPtc[2] = tpsaPtcTemp[2];
        tpsaPtc[3] = tpsaPtcTemp[3];
    }
    else if (xk < 0.0) {
        sqk = sqrt(-outin);
        co = cos(sqk * tau);
        si = sin(sqk * tau);
        ch = cosh(sqk * tau);
        sh = sinh(sqk * tau);

        tpsaPtcTemp[0] = tpsaPtc[0] * ch + 
                         tpsaPtc[1] * sh / sqk;
        tpsaPtcTemp[1] = tpsaPtc[0] * sh * sqk + 
                         tpsaPtc[1] * ch;
        tpsaPtcTemp[2] = tpsaPtc[2] * co + 
                         tpsaPtc[3] * si / sqk;
        tpsaPtcTemp[3] = -tpsaPtc[2] * si * sqk + 
                         tpsaPtc[3] * co;

        tpsaPtc[0] = tpsaPtcTemp[0];
        tpsaPtc[1] = tpsaPtcTemp[1];
        tpsaPtc[2] = tpsaPtcTemp[2];
        tpsaPtc[3] = tpsaPtcTemp[3];
    }
}

