// CODE TO Calculate bond yield via newton's method


#include <stdio.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double bond_price (double v_cash_flow[], double t_cash_flow[], double x_old, int N)

{
    double B=0;
    for(int i = 0;i < N;i++){
        B = B +(v_cash_flow[i]*exp(-1*t_cash_flow[i]*x_old));}
    return B;
}


double bond_derivative_price (double v_cash_flow[], double t_cash_flow[], double x_old, int N)
{
    double B_D = 0;
    for(int i = 0;i < N;i++){
        B_D = B_D + t_cash_flow[i] * v_cash_flow[i] * exp(-1 * x_old * t_cash_flow[i]);}
    return B_D;
}



int main() {
    int N = 7;
    // ENTER Number of cash flows
    
    double t_cash_flow[] = {0.33333333333333, 0.83333333333333, 1.33333333333333, 1.83333333333333, 2.33333333333333, 2.83333333333333, 3.33333333333333};
    // array of cash flow dates MAKE SURE THESE ARE CORRECT AND ASSCENDING
    
    double v_cash_flow[] = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 104.5};
    //array of cash flows DOUBLE CHECK THESE
    
    double price = 113.4551703428;
    // enter price here
    
    double x_0 = 0.1; // initial guess
    
    double x_new = x_0;
    double x_old = x_0 - 1;
    
    double tol = pow(10, -6);
    
    
    
    while (abs(x_new - x_old) > tol)
        
    {
        cout << std::setprecision(6)  << x_new << endl;
        x_old = x_new;
        x_new = x_old + ((bond_price(v_cash_flow, t_cash_flow, x_old, N)-price) / bond_derivative_price(v_cash_flow, t_cash_flow, x_old, N));
        
    }
    
    
    
    cout << "Yield is: " << std::setprecision(6)  << x_new << endl;
    
    
    return 0;
}
