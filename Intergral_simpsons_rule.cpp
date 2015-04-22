//
//  main.cpp
//  Midpoint rule - intergrals - for FINAL
//
//  Created by Gabriel Levine on 11/27/14.
//  Copyright (c) 2014 Gabriel Levine. All rights reserved.
//

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
# include <tgmath.h>
#include <iomanip>
using namespace std;

// insert function intergral here
double func1(double x){
    return exp(-1*pow(x, 2));
    
}
/*
double midpoint(double a, double b, int n){
    int i;
    double h = ((b - a)/n);
    double I_midpoint = 0;
    for(i = 1;i <= n;i++)
    {
        I_midpoint = I_midpoint + func1(a+((i - 0.5)*h));
        
    }
    I_midpoint = h * I_midpoint;
    
    cout << "interval: " << n << ", "<< std::setprecision(8) << I_midpoint << endl;
    return I_midpoint;
}
*/

double simpsons(double a, double b, int n){
    int i;
    double h = ((b - a)/n);
    double I_simpson = (func1(a)/6) + (func1(b)/6);
    for(i = 1;i <= n-1;i++)
    {
        I_simpson = I_simpson + func1(a + i*h)/3;
        
    }

    for(i = 1;i <= n;i++)
    {
        I_simpson = I_simpson + 2*func1(a+((i - 0.5)*h))/3;
        
    }
    
    I_simpson = h * I_simpson;
    
    cout << "interval: " << n << ", "<< std::setprecision(8) << I_simpson << endl;
    return I_simpson;
}


int main() {
    double a = 0; //left endpoint
    double b = 2; // right endpoint
    
    double tol = 5*pow(10, -7); // set tol
    
    
    int n = 4; // number of intervals, initial
    
    double I_old = simpsons(a,b,n);
    n= 2*n;
    double I_new = simpsons(a, b, n);
    while (fabs(I_new - I_old) > tol)
    {
        I_old = I_new;
        n= 2*n;
        I_new = simpsons(a, b, n);
    }
    double I_approx = I_new;
    
    cout << "final interval: " << n << ", "<< std::setprecision(8) << I_approx << endl;
    
    
    
    return 0;
}
