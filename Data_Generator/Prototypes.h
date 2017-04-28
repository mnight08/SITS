#include <complex>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
#ifndef Prototypes
#define Prototypes
class point{
	public:
	double x, y;
	inline void set(double xp,double yp){x=xp;y=yp;}
	point(){x=0;y=0;}
	point(double xp, double yp){x=xp;y=yp;}
};

//!!!!!!!Functions to data Generation!!!!!!!!!!!!!!!!!!!!!//
inline int Field_On_Support_Index(int s, int w, int j);
void Read_In_Reflectivity(double *reflectivity,string file);
int Size_Of_Support(double *reflectivity);
complex<double> Incident_Field(int s, int w, double x, double y);
double Norm(double x,double y, double z);
inline complex<double> G(double x,double y, double z,double w);
complex<double> Scattered_Field_At_Point(double x, double y,int w,int s,double* reflectivity,complex<double>* field_on_target_support,int S);
complex<double> *Generate_Data(double *reflectivity,int K);
inline void Output_Data(complex<double> *data, string file_name);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//




//!!!!!!!!!!!!!!!!!!Functions for image reconstruction!!!!!!//
inline int Grid_Index( int Box,int t, int l);
inline point Grid_Point(double * grids, int Box, int t1, int t2, int l);
inline void Proto_Grid(double *fillme);
inline void Grid(double *protogrid, int Box, int l);
inline void Build_Grids(double *&grids);
inline double L1d(double *grids, int l, int Box, int t,double x);
inline double L2d(double *grids, int t1, int t2, int B, int l,point y);
inline void Swap(complex<double>*& current, complex<double>* &previous);
inline void Zero(complex<double> * weights);
inline double Center1d(int l, int box);
inline point Center(int l, int box);
inline int A(int n, int l);
inline int B(int n, int l);
inline int Parent(int A);
inline int Child(int B,int c);
inline int Weight_Index(int A, int B, int l, int t1, int t2);
inline int Box(point p, int l);
inline int Image_Index(int m, int n);
inline int Data_Index(int s,int w);
inline double X(int m);
inline double Y(int n);
inline double Gamma_D_X(double sl);
inline double Gamma_D_Y(double sl);
//normalized frequency
inline double Q(int e);
inline double W(int w);
inline double S(int s);
void Image_Recovery_Butterfly(double *& recovered_reflectivity, complex<double>* data);
void Image_Recovery_Direct(double * recovered_reflectivity, complex<double>* data);
void Output_Reflectivity(double * reflectivity,string file_name);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

//!!!!!!!!!!!!!!!!!!!Input functions!!!!!!!!!!!!!!!!!!!!!!!!!!!//
complex<double> f(point y,complex<double>*d);
complex<double> p(double w);
complex<double> interp(complex<double>* d, double sl,double qj);
double phi(point x,point y);
double Gamma_X(double s);
double Gamma_Y(double s);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
inline void Test_Grid_Creation();
inline void Test_Indexing();
inline void Test_Image_Loading();
	//Might not be able to do these here.
inline void Test_Data_Generation();
inline void Test_Data_Loading();
inline void Test_Image_Recovery();
inline void Test_Image_Output();

inline void Print_Grids(double *grids);
inline void Print_Weights(complex<double>* weights,int l);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
#endif
