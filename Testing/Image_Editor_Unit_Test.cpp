#include "input.h"
#include "data_generation.h"
#include "image_recovery.h"
#include <ofstream>
double phi(double x1,double x2, double y1,double y2);
double f(double y1, double y2);
double psi(double x1,double x2, double y1,double y2);

//This file will test each function separately to make sure they are
//working correctly. Testing data will be output to a file.
//The printing function live here since they are strictly for testing purposes



int main(int argc, char *argv[])
{

	ofstream log;
	log.open("log.txt");
	Test_Grid_Creation(log);
	Test_Chebyshev_Interpolation(ostream);
	Test_Indexing(ostream);
	Test_Image_Loading(ostream);
	//Might not be able to do these here.
	Test_Data_Generation(ostream);
	Test_Data_Loading(ostream);
//	Test_Image_Recovery();
///	Test_Image_Output();
	


//	Read_In_Reflectivity(reflectivity,"EasyLargest.in");
//	Output_Reflectivity(reflectivity,"test.out");
//	cout<<"Reading in Image"<<endl;
//	Read_In_Reflectivity(reflectivity,"EasyLargest.in");
//	Output_Reflectivity(reflectivity,"test.out");
//	cout<<"Reading in data"<<endl;
//	Read_In_Data(data,"data.in");
//	cout<<"Read in data"<<endl;
	//Reconstruct image.  Need to specify phase, and amplitude functions in input file
	//Initialize the recovered image
//	double *recovered_reflectivity=new double[M];

//	cout<<Size_Of_Support(reflectivity);
//	std::cout<<"Generating Data"<<endl;
//	data=Generate_Data(reflectivity,K);
//	cout<<"Data was generated"<<endl;
//	Output_Data(data,"dataEasyLargest.in");

	std::cout<<"Testing is complete."<<endl;

	return 0;
}


inline void Output_Weights(ostream os,complex<double>* weights,int l)
{
	for(int n=0;n<M;n++){
			//go through the children of B paired with the parent of A to compute weights
			for(int t1=0;t1<q;t1++){
				for(int t2=0;t2<q;t2++){
					os<<weights[Weight_Index(A(n,l),B(n,l),l,t1,t2)]<<" ";
				}
			}
	}
}

inline void Output_Grids(ostream os,double *grids)
{
	for( int l=0;l<=L;l++){
		//go through each offset at that scale.
		for(int Box=0;Box<(int)(pow(2.0,l));Box++){
			for( int t=0;t<q;t++)
			{
				os<<grids[Grid_Index(Box,t,l)]<<" ";
			}
			os<<endl;
		}
		os<<endl;
	}
	cout<<endl<<endl;

}

inline void Test_Grid_Creation(){
	double *grids;
	Build_Grids(grids);
	Ouput_Grids(log,grids);
}



inline void Test_Indexing(ofstream log){
	
	double *indices;
	Build_Grids(indices);
	Ouput_Grids(log,indices);
}
inline void Test_Image_Loading(){

}
	//Might not be able to do these here.
inline void Test_Data_Generation(){

}
inline void Test_Data_Loading(){

}
inline void Test_Image_Recovery(){

}
inline void Test_Image_Output(){
	
}


inline void Test_Chebyshev_Interpolation(){
	
}
