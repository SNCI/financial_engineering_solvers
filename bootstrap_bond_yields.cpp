// CODE TO Calculate bond yield via newton's method -- BOOTSTRAPPING!!!!
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
// r(0, 1.5) = 0.023752

// enter non-linear function set to 0 here.
double Func1(double x)
{
    double temp = -116 + 4*exp(-0.5*0.017245758276) + 4*exp(-0.019491516552) +4*exp(-1.5*0.023752) + 4*exp(-2*(0.2*x +0.8*0.023752)) + 4*exp(-2.5*(0.4*x + 0.6*0.023752)) +4*exp(-3*(0.6*x + 0.4*0.023752)) + 4*exp(-3.5*(0.8*x + 0.2*0.023752)) + 104*exp(-4*x)  ;
    
    return temp;
}

// enter derivative of that funtion here
double DerivativeFunc1(double x)
{
    double d_temp = 4*-2.5*0.4*exp(-2.5*(0.4*x + 0.6*0.023752)) + 4*-3*0.6*exp(-3*(0.6*x + 0.4*0.023752)) + 4*-3.5*0.8*exp(-3.5*(0.8*x + 0.2*0.023752)) + 104*-4*exp(-4*x);
    
    return d_temp;
}



int main() {
    
    
    
    double x_0 = 0.05;
    // initial guess
    
    double x_new = x_0;
    double x_old = x_0 - 1;
    
    double tol = pow(10, -6);
    // to change tol
    
    while (abs(x_new - x_old) > tol)
        
    {
        cout << std::setprecision(6)  << x_new << endl;
        x_old = x_new;
        x_new = x_old - (Func1(x_old)/DerivativeFunc1(x_old));
        
    }
    
    
    
    cout << "Solution is: " << std::setprecision(6)  << x_new << endl;
    cout << " " << endl;
    double rate1 = (0.2)*(0.0373664) + (0.8)*(0.023752); cout << "2: "<< std::setprecision(6) << rate1 << endl;
    double rate2 = (0.4)*(0.0373664) + (0.6)*(0.023752); cout << "2.5: "<< std::setprecision(6) << rate2 << endl;
    double rate3 = (0.6)*(0.0373664) + (0.4)*(0.023752); cout << "3: "<< std::setprecision(6) << rate3 << endl;
    double rate4 = (0.8)*(0.0373664) + (0.2)*(0.023752); cout << "3.5: "<< std::setprecision(6) << rate4 << endl;


    
    
    return 0;
}
