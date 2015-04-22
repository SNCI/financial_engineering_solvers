//
//  CODE to find implied Volatility

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

double d_1(double t, double S, double K, double T, double v, double q, double r) {
    return (log(S/K) + ((r - q + pow(v,2)/2)*T))/(v*sqrt(T));
    
}

double d_2(double t, double S, double K, double T, double v, double q, double r) {
    return (log(S/K) + ((r - q - pow(v,2)/2)*T))/(v*sqrt(T));
}


double put_price_norm(double t, double S, double K, double T, double v, double q, double r) {
    double N_d1 = norm_cdf(-1*d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(-1*d_2(t, S, K, T, v, q, r));
    return (K * exp(-r*T)*(N_d2)) - (S * exp(-q*T) * (N_d1));
}


double call_price_norm(double t, double S, double K, double T, double v, double q, double r) {
    double N_d1 = norm_cdf(d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(d_2(t, S, K, T, v, q, r));
    return (S * (exp(-q*T) * (N_d1)) - (K * exp(-r*T) * N_d2));
}

double call_price(double t, double S, double K, double T, double v, double q, double r) {
    double N_d1 = SimpsonsDistribution(d_1(t, S, K, T, v, q, r));
    double N_d2 = SimpsonsDistribution(d_2(t, S, K, T, v, q, r));
    return (S * (exp(-q*T) * (N_d1))) - (K * exp(-r*T) * N_d2);
}



double put_price(double t, double S, double K, double T, double v, double q, double r) {
    double N_d1 = SimpsonsDistribution(-1*d_1(t, S, K, T, v, q, r));
    double N_d2 = SimpsonsDistribution(-1*d_2(t, S, K, T, v, q, r));
    return (K * exp(-r*T)*(N_d2)) - (S * exp(-q*T) * (N_d1));
}

double delta_call(double t, double S, double K, double T, double v, double q, double r)
{
    double d1 = (log(S/K) + (r - q + ((v*v)/2))*(T-t))/(v*sqrt(T-t));
    double delta_c =(exp(-1*q*T-t) * norm_cdf(d1));
    return delta_c;
}

double delta_put(double t, double S, double K, double T, double v, double q, double r)
{
    double d1 = (log(S/K) + (r - q + ((v*v)/2))*(T-t))/(v*sqrt(T-t));
    double delta_c =(-1*exp(-1*q*T-t) * norm_cdf(-1*d1));
    return delta_c;
}

double gamma_call(double t, double S, double K, double T, double v, double q, double r)
{
    double d1 = (log(S/K) + (r - q + ((v*v)/2))*(T-t))/(v*sqrt(T-t));
    double gamma_c = (exp(-1*q*(T-t))/(K*v*sqrt(T-t))) * ((1.0/sqrt(2*M_PI))) * (exp((d1*d1)/2));
    return gamma_c;
}

double vega_call(double t, double S, double K, double T, double v, double q, double r)
{
    return S*exp(-q*(T-t))*sqrt(T-t)*exp(-1*d_2(t, S, K, T, v, q, r)*0.5)/sqrt(2*M_PI);
}

double theta_call(double t, double S, double K, double T, double v, double q, double r)
{
    double N_d1 = norm_cdf(d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(d_2(t, S, K, T, v, q, r));
    double temp = (-1*(S*v*exp(-q*(T-t)))/(2*sqrt(2*M_PI*(T-t)))*exp(-1*d_1(t, S, K, T, v, q, r)*d_1(t, S, K, T, v, q, r)/2)) + (q*S*exp(-q*(T-t))*N_d1) - (r*K*N_d2*exp(-r*(T-t))) ;
    return temp;
}

double theta_put(double t, double S, double K, double T, double v, double q, double r)
{
    double N_d1 = norm_cdf(-1*d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(-1*d_2(t, S, K, T, v, q, r));
    double temp = (-1*(S*v*exp(-q*(T-t)))/(2*sqrt(2*M_PI*(T-t)))*exp(-1*d_1(t, S, K, T, v, q, r)*d_1(t, S, K, T, v, q, r)/2)) - (q*S*exp(-q*(T-t))*N_d1) + (r*K*N_d2*exp(-r*(T-t))) ;
    return temp;
}

double rho_call(double t, double S, double K, double T, double v, double q, double r)
{
    //   double N_d1 = SimpsonsDistribution(d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(d_2(t, S, K, T, v, q, r));
    double temp = N_d2*K*(T-t)*exp(-r*(T-t));
    return temp;
}

double rho_put(double t, double S, double K, double T, double v, double q, double r)
{
    //   double N_d1 = SimpsonsDistribution(d_1(t, S, K, T, v, q, r));
    double N_d2 = norm_cdf(-1*d_2(t, S, K, T, v, q, r));
    double temp = N_d2*-1*K*(T-t)*exp(-r*(T-t));
    return temp;
}

double DerivativeFunc1(double t, double S, double K, double T, double v, double q, double r)
{
    // Just the derivative, need for Newton's method
    double N_d2 = norm_cdf(d_2(t, S, K, T, v, q, r));
    double d_temp = 1 - N_d2*exp(-r*(T-t));
    
    return d_temp;
}
double d_call_interms_K(double t, double S, double K_old, double T, double v, double q, double r) {
    //  double N_d1 = SimpsonsDistribution(d_1(t, S, K_old, T, x, q, r));
    double N_d2 = norm_cdf(d_2(t, S, K_old, T, v, q, r));
    return -1*exp(-r*T) * N_d2;
}




int main() {
    // First we create the parameter list
    double t = 0.0; // start time
    double T = 0.5; // End time of option in terms of Year 1.0
    
    double K = 40.0; // strike
    double S = 42.0; // spot price
    
    double v = 0.3; // sigma, volitility
    
    double q = 0.03; // dividends
    double r = 0.05; // interest rate
    
    
    double call = call_price(t, S, K, T, v, q, r); // make sure you have the right function, norm or simpsons!!!
    double put = put_price(t, S, K, T, v, q, r);
    double d1_norm= norm_cdf(d_1(t, S, K, T, v, q, r));
    double d2_norm= norm_cdf(d_2(t, S, K, T, v, q, r));
    double d1_simp= SimpsonsDistribution(d_1(t, S, K, T, v, q, r));
    double d2_simp= SimpsonsDistribution(d_2(t, S, K, T, v, q, r));
    double delta_c = delta_call(t, S, K, T, v, q, r);
    double delta_p = delta_put(t, S, K, T, v, q, r);
    double gamma_cp = gamma_call(t, S, K, T, v, q, r);
    double vega_cp = vega_call(t, S, K, T, v, q, r);
    double theta_c = theta_call(t, S, K, T, v, q, r);
    double theta_p = theta_put(t, S, K, T, v, q, r);
    double rho_c = rho_call(t, S, K, T, v, q, r);
    double rho_p = rho_put(t, S, K, T, v, q, r);
    
    cout << "Call Price: " << std::setprecision(11)  << call << endl;
    cout << "Put Price: " << std::setprecision(11)  << put << endl;
    cout << "  " << endl;
    
    cout << "N(d1) cum_dist_norm: " << std::setprecision(11)  << d1_norm << endl;
    cout << "N(d2) cum_dist_norm: " << std::setprecision(11)  << d2_norm << endl;
    cout << "N(d1) simpsons: " << std::setprecision(11)  << d1_simp << endl;
    cout << "N(d2) simpsons: " << std::setprecision(11)  << d2_simp << endl;
    cout << "  " << endl;
    
    cout << "Delta Call: " << std::setprecision(11)  << delta_c << endl;
    cout << "Delta Put: " << std::setprecision(11)  << delta_p << endl;
    cout << "Gamma Call/Put: " << std::setprecision(11)  << gamma_cp << endl;
    cout << "Vega Call/Put: " << std::setprecision(11)  << vega_cp << endl;
    cout << "Theta Call: " << std::setprecision(11)  << theta_c << endl;
    cout << "Theta Put: " << std::setprecision(11)  << theta_p << endl;
    cout << "Rho Call: " << std::setprecision(11)  << rho_c << endl;
    cout << "Rho Put: " << std::setprecision(11)  << rho_p << endl;
    
    
    
    
    
    return 0;
}



