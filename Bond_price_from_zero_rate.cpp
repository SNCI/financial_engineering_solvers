#include <stdio.h>
#include <stdlib.h>
# define _USE_MATH_DEFINES
#include <math.h>
# include <iostream>
# include <ostream>
# include <iomanip>
using namespace std;


// CODE to compute the PRICE and discounts -- given ZERO RATE CURVE


int main()

{

    int n=7; // ENTER number of cash flows
    
    int i;
    
    double B=0; // price of the bond
    
    double t_cash_flow[] = {0.33333333333333, 0.83333333333333, 1.33333333333333, 1.83333333333333, 2.33333333333333, 2.83333333333333, 3.33333333333333};
    // array of cash flow dates MAKE SURE THESE ARE CORRECT AND ASSCENDING
    
    double v_cash_flow[] = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 104.5};
    //array of cash flows DOUBLE CHECK THESE
    
    for(i = 0;i < n;i++)
    {   double r_zero = 0.04 + ((log(1+2*t_cash_flow[i]))/200);
        // INSERT NEW zero rate above. Double check!!!!
        
        double disc = exp(-1* t_cash_flow[i] * r_zero); // stays the same
        cout << "Disc: " << std::setprecision(6) << disc << ", time " << t_cash_flow[i] << endl;
        B = B + v_cash_flow[i] * disc;
        
    }
    cout << "Bond Value is: " << std::setprecision(13) << B << endl;
    
    
    
    return 0;
}
