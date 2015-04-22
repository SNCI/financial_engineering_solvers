// CODE TO Calculate bond yield via newton's method -- BOOTSTRAPPING!!!!
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


// enter non-linear function set to 0 here.
double Func1(double x)
{
    double temp = -103 + 2.5*exp(-0.5*(0.5*x + 0.0075)) +102.5*exp(-x);
    
    return temp;
}

// enter derivative of that funtion here
double DerivativeFunc1(double x)
{
    double d_temp = 2.5*(-0.25)*exp(-0.5*(0.5*x + 0.0075)) - 102.5*exp(-x);
    
    return d_temp;
}



int main() {
   

    
    double x_0 = 0.002;
    // initial guess
    
    double x_new = x_0;
    double x_old = x_0 - 1;
    
    double tol = pow(10, -6);
    // to change tol
    
    while (abs(x_new - x_old) > tol)
        
    {
        cout << std::setprecision(9)  << x_new << endl;
        x_old = x_new;
        x_new = x_old - (Func1(x_old)/DerivativeFunc1(x_old));
        
    }
    
    
    
    cout << "Solution is: " << std::setprecision(11)  << x_new << endl;
    cout << " " << endl;
    double rate1 = (0.5)*(0.019491516552) + (0.5)*(0.015); cout << "1: "<< rate1 << endl;
    // double rate2  = (0.33333333333333)*(0.024964739601) + (0.66666666666667)*(0.02381303); cout <<"2: "<< rate2 << endl;
    
    return 0;
}
