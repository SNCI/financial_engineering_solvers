#include <stdio.h>
#include <stdlib.h>
# define _USE_MATH_DEFINES
#include <math.h>
# include <iostream>
# include <ostream>
# include <iomanip>
// # include "pointclass.h"
using namespace std;
// CODE to compute the PRICE, DURATION and CONVEXITY of a bond given the yield
int main()
{
    int i;
    
    int n=7; // number of cash flows
    
    double B=0; // price of the bond
    double D=0;
    double C=0;
    double dollar_D=0;
    double dollar_C=0;
    double dv01=0;
   

    double disc[n];
    
    double B_price = 113.4551703428; // copy from other output

    
    double y = 0.0499095 ; // yeild of the bond COPY FROM output
    
    
    double t_cash_flow[] = {0.33333333333333, 0.83333333333333, 1.33333333333333, 1.83333333333333, 2.33333333333333, 2.83333333333333, 3.33333333333333};
    // array of cash flow dates MAKE SURE THESE ARE CORRECT AND ASSCENDING
    
    double v_cash_flow[] = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 104.5};
    //array of cash flows DOUBLE CHECK THESE
    
    for(i = 0;i < n;i++)
    {
        disc[i] = exp(-1*t_cash_flow[i]*y);
        B = B +(v_cash_flow[i]*disc[i]);
        D = D + (t_cash_flow[i]*v_cash_flow[i]*disc[i]);
        C = C + ((t_cash_flow[i]*t_cash_flow[i])*v_cash_flow[i]*disc[i]);
    }
    D = D/B;
    C = C/B;
    
    dollar_D = B_price*D;
    dollar_C = B_price*C;
    dv01 = dollar_D/10000;

    
    cout << "for a Bond with yeild " << y << endl;
//    cout << "Bond Value is: " << std::setprecision(10) << B << endl;
    cout << "Modified Duration is: " << std::setprecision(7) << D << endl;
    cout << "Dollar Duration is: " << std::setprecision(10) << dollar_D << endl;
    cout << "Bond Convexity is: " << std::setprecision(10) << C << endl;
    cout << "Dollar Convexity is: " << std::setprecision(11) << dollar_C << endl;
    cout << "DV01 is: " << std::setprecision(10) << dv01 << endl;

    
    return 0;
}
