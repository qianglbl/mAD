//
// Created by Ji Qiang.
//

#include "TPSAad.h"
#include <iostream>

int main(){
    TPSAad x1,x2,x3,x4,x5,func;

    x1.assign(3.0,1);
    x2.assign(0.1,2);
    x3.assign(0.1,3);
    x4.assign(1.0,4);
    x5.assign(1.0,5);

    //test function
    func = (2*cos(x1/x2)+x1/x2+exp(x2))*x3+2*x4+sinh(x5);

    cout<<"func. value and its derivatives w.r.t. 5 variables"<<endl;
    cout<<func.map[0]<<" "<<func.map[1]<<" "<<func.map[2]<<" "<<func.map[3]<<" "<<func.map[4]<<" "<<func.map[5]<<endl;
}

