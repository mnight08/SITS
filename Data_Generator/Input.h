#include "prototypes.h"
#ifndef INPUT
#define INPUT
using namespace std;
//Parameters
const complex<double> i(0,1);
const double pi=3.14159265359;
const double c_0=299792458;
//Assume the samples are uniform and given at $(i+1/2N,j+1/2N)$ for $i,j=0,1,\ldots N-1$
//N is the number of slow time, frequency samples.  There are $N^2$ space samples
const int N=4;
const int M=N*N;
//number of scattering events to consider.
const int K=1;
//center frequency
const double w_0=pow(10,9);
//Bandwidth
const double Band=pow(10,8);
//radius of flight path
const double r=100;
//Flight height
const double H=100;
//Be careful with this since, you assumed L=Log2(N)
const int L=log2(N);
const int q=16;
const int r_eps=q*q;
//image plane height
const double image_height=0;
//Number of weights at each level.
const int P=N*N*q*q;
const string filename="Easy";

inline complex<double> G(double x,double y, double z,double w){
	return exp(-i*w*(Norm(x,y,z))/c_0)/(4*pi*Norm(x,y,z));
}

inline double Norm(double x,double y, double z){
	return sqrt(x*x+y*y+z*z);
}

double Gamma_X(double s){
	return r*cos(2*pi*s);
}

double Gamma_Y(double s){
	return r*sin(2*pi*s);
}

//Phase function.  All arguments are in $[0,1]^2$
double phi(point x,point y){
	return -(w_0-Band/2+y.y*Band)*(2*(Gamma_X(y.x)*x.x+Gamma_Y(y.x)*x.y-r)/c_0);
}


//Linear interpolation.  Assume that data is sorted in frequency and slow time.
complex<double> interp(complex<double>* d, double sl,double qj){	
	int s=(int)(sl*N);
	int w=(int)(qj*N);
	return d[Data_Index(s,w)];
}

complex<double> p(double w){
	if((w<w_0+Band/2)&&(w>w_0-Band/2))
		return 1;
	else 
		return 0;
}
	

//normalized function domain
complex<double> f(point y,complex<double>*d){
	return (conj(p(w_0-Band/2+y.y*Band))*interp(d,y.x,y.y))
		/((w_0-Band/2+y.y*Band)*abs(p(w_0-Band/2+y.y*Band))*abs(p(w_0-Band/2+y.y*Band)));
}
#endif
