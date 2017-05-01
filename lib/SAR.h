
double Gamma_X(double s){
	return r*cos(2*pi*s);
}

double Gamma_Y(double s){
	return r*sin(2*pi*s);
}

//Phase function.  All arguments are in $[0,1]^2$
double phi(point x,point y){
	return -(w_0-Band/2+y.y*Band)*(2*(Gamma_X(y.x)*x.x+Gamma_Y(y.x)*x.y-r)/c_0);
inline complex<double> G(double x,double y, double z,double w){
	return exp(-i*w*(Norm(x,y,z))/c_0)/(4*pi*Norm(x,y,z));
}


complex<double> p(double w){
	if((w<w_0+Band/2)&&(w>w_0-Band/2))
		return 1;
	else 
		return 0;
}
	

