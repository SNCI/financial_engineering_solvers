//
//  CODE to find implied Volatility for put option.

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// Standard normal probability density function
double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// An approximation to the cumulative distribution function for the standard normal distribution
long double norm_cdf(double t) {
    double nn;
    double z = fabs(t);
    double et = exp(-1*(t*t/2));
    double y = 1.0/(1.0 + 0.2316419*z);
    double a1 = 0.319381530;
    double a2 = -0.356563782;
    double a3 = 1.781477937;
    double a4 = -1.821255978;
    double a5 = 1.330274429;
    double m = 1.0 - (et * ((a1 * y + a2 * pow(y,2) + a3 * pow(y,3) + a4 * pow(y,4) + a5 * pow(y,5))) / (sqrt(2*M_PI)));
    
    if (t >= 0.0) {
        return nn = m;
    } else {
        return nn = 1 - m;
    }
}

//  This is the Statistical ERROR funtion
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}


// calculate the value of the cumulative normal distribution using Simpson's method
double SimpsonsDistribution(double x)
{
    if(x<-10.)return 0.; // return sensible values on limits
    if(x>10.)return 1.;
    // number of steps
    int N=1000;
    // range of integration
    double a=0,b=x;
    // local variables
    double s,h,sum=0.;
    // initialise the variables
    h=(b-a)/N;
    // add in the first few terms
    sum = sum + exp(-a*a/2.) + 4.*exp(-(a+h)*(a+h)/2.);
    // and the last one
    sum = sum + exp(-b*b/2.);
    // loop over terms 2 up to N-1
    for(int i=1;i<N/2;i++)
    {
        s = a + 2*i*h;
        sum = sum + 2.*exp(-s*s/2.);
        s = s + h;
        sum = sum + 4.*exp(-s*s/2.);
    }
    // complete the integral
    sum = 0.5 + h*sum/3./sqrt(8.*atan(1.));
    // return result
    return sum;
}

double d_1(double t, double S, double K, double T, double x, double q, double r) {
    return (log(S/K) + ((r - q + pow(x,2)/2)*T))/(x*sqrt(T));
    
}

double d_2(double t, double S, double K, double T, double x, double q, double r) {
    return (log(S/K) + ((r - q - pow(x,2)/2)*T))/(x*sqrt(T));
}





double call_price_norm(double t, double S, double K, double T, double v, double q, double r) {
    double d1_cdf;
    double d2_cdf;
    double time = T-t;
    d1_cdf = phi(d_1(t, S, K, T, v, q, r));
    d2_cdf = phi(d_2(t, S, K, T, v, q, r));
    return (S * (exp(-q*time) * (d1_cdf)) - (K * exp(-r*time) * d2_cdf));
}

double call_price_norm_2(const double& t, const double& S, const double& K, const double& T, const double& v, const double& q, const double& r) {
    double d1_cdf;
    double d2_cdf;
    double time = T-t;
    d1_cdf = norm_cdf(d_1(t, S, K, T, v, q, r));
    d2_cdf = norm_cdf(d_2(t, S, K, T, v, q, r));
    return (S * (exp(-q*time) * (d1_cdf)) - (K * exp(-r*time) * d2_cdf));
}

double put_price_norm(const double& t, const double& S, const double& K, const double& T, const double& x, const double& q, const double& r) {
    double N_d1 = norm_cdf(-1*d_1(t, S, K, T, x, q, r));
    double N_d2 = norm_cdf(-1*d_2(t, S, K, T, x, q, r));
    return (K * exp(-r*T)*(N_d2)) - (S * exp(-q*T) * (N_d1));
}


double call_price(double t, double S, double K, double T, double x, double q, double r) {
    double d1_cdf;
    double d2_cdf;
    d1_cdf = SimpsonsDistribution(d_1(t, S, K, T, x, q, r));
    d2_cdf = SimpsonsDistribution(d_2(t, S, K, T, x, q, r));
    return (S * (exp(-q*T) * (d1_cdf)) - (K * exp(-r*T) * d2_cdf));
}

double put_price(double t, double S, double K, double T, double x, double q, double r) {
    double N_d1 = SimpsonsDistribution(-1*d_1(t, S, K, T, x, q, r));
    double N_d2 = SimpsonsDistribution(-1*d_2(t, S, K, T, x, q, r));
    return (K * exp(-r*T)*(N_d2)) - (S * exp(-q*T) * (N_d1));
}

double vega_call(double t, double S, double K, double T, double x, double q, double r)
{
    double d1 = (log(S/K) + (r - q + ((x*x)/2))*(T))/(x*sqrt(T-t));
    return (1.0/sqrt(2*M_PI))*sqrt(T)*(exp(-1*(pow(d1, 2)/2)))*(S*exp(-q*T));
}



int main() {
    double t = 0.0; // start time
    double T = 0.58333333333333; // End time of option in terms of Year 1.0
    double K = 52; // strike price
    double P = 5.06; // put price
    double q = 0.0075; // dividends
    double S = 50.0; // spot price
    double r = 0.0175; // interest rate
    
    
    double x_0 = 0.2; // sigma, volitility initial guess

    double x_new = x_0;
    double x_old = x_0 - 1;
    
    double tol = pow(10, -8);

    
    while (fabs(x_new - x_old) > tol)
        
    {
        cout << std::setprecision(8)  << x_new << endl;
        x_old = x_new;
        x_new = x_new - ((put_price_norm(t,S,K,T,x_new,q,r) - P) / vega_call(t,S,K,T,x_new,q,r));
    }
    
    cout << "Implied volatility is: " << std::setprecision(8)  << x_new << endl;
    
    //   double imp_approx = ((sqrt(2*M_PI)) / (S*sqrt(T))) * ((C-(((r-q)*T)/2)*S) / ((1-((r+q)*T*0.5))));
    // double error = (fabs(imp_approx - x_n)) / x_n;
    //cout << "Relative Error: " << std::setprecision(13)  << error << endl;
    
    
    return 0;
}



