#include "input.h"
//Input  and output````````````````````````````///
double* Read_In_Reflectivity(string infile_name){

    
    //Reflectivities are real.  There are N^2=M image points.
    double *reflectivity=new double[M];
    ifstream input;
    input.open((infile_name+".in").c_str());
    for(int n=0;n<N;n++) {
	for(int m=0;m<N;m++){
	   input>>reflectivity[Image_Index(m,n)];
	}
    }
    input.close();
}

inline void Output_Data(complex<double> *data, string outfile_name){
	ofstream output;
	output.open((outfile_name+".out").c_str());
	for(int s=0;s<N;s++){
		for(int w=0;w<N;w++){
		   output<<data[Data_Index(s,w)]<<" ";
		}
		output<<endl;
	}
	output.close();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~///

//Number of non zero entries in $v$
int Size_Of_Support(double *reflectivity){
	int size=0;
	for(int m=0;m<N;m++)
		for(int n=0;n<N;n++)
			if(reflectivity[Image_Index(m,n)]!=0)
				size++;
	return size;
}

//Incident field at a point, slow time sample, and frequency sample. Assumes that Antenna height is fixed.
complex<double> Incident_Field(int s, int w, double x, double y,double z){
	return G(x-Gamma_X(S(s)),y-Gamma_Y(S(s)),z-H,W(w));
}

///Evaluate scattering operator at a point
complex<double> Scattered_Field_At_Point(double x, double y,int w,int s,double* 
		reflectivity,complex<double>* field_on_target_support,int S){
	complex<double> field_at_point(0,0);
	//index for non zero points of reflectivity
	int j=0;
	//integrate over target support
	for(int m=0;m<N;m++){
		for(int n=0;n<N;n++){
			if(reflectivity[Image_Index(m,n)]!=0){
					if(X(m)!=x||Y(n)!=y){
						field_at_point+=field_on_target_support[Field_On_Support_Index(s,w,j)]*G(x-X(m),y-Y(n),0,W(w));
					}
				j++;
			}
		}
	}
	return (-1.0*W(w)*W(w))*field_at_point/((double) N*N);
}
 
void Generate_Data(string infilename){
    double *reflectivity=Read_In_Reflectivity(infilename);
    complex<double> *data=Generate_Data(reflectivity,1);
    string outfilename=infilename.erase(infilename.end()-3,infilename.end())+".out";
    Output_Data(data,outfilename);
}

//Generate data due to the reflectivity funciton and number of scattering events K
complex<double> *Generate_Data(double *reflectivity,int K)
{
    //Data is complex valued.  There are N^2=M data points.
    //complex<double> *data=new complex<double>[M];
    
	complex<double> *data=new complex<double>[N*N];
  		
	int SP=Size_Of_Support(reflectivity);
	std::cout<< "Size of Support is "<<SP<<endl;
	
	//Initialize to the incident field.
	complex<double>	*field_on_target_support=new complex<double>[N*N*SP];

	//target support index
	int j=0;
	//Generate incident field on target support
	for(int n=0;n<N;n++){
		for(int m=0;m<N;m++){
			if(reflectivity[Image_Index(m,n)]!=0){
				//Generate field on supp for each pulse
				for(int s=0;s<N;s++){
					//Generate field on supp for each frequency
					for(int w=0;w<N;w++){
						field_on_target_support[Field_On_Support_Index(s,w,j)]=Incident_Field(s,w,X(m),Y(n),0);
						}
				}
				j++;
			}
		}
	}
	std::cout<<"Generated incident field on target support"<<endl;
	for (int k=1;k<K;k++) {
		//Scatter previous field on target support
		j=0;
		std::cout<<"You are in the K loop"<<endl;
		for(int n=0;n<N;n++){
			for(int m=0;m<N;m++){
				if(reflectivity[Image_Index(m,n)]!=0){
					//Generate field on supp for each pulse
					for(int s=0;s<N;s++){
						//Generate field on supp for each frequency
						for(int w=0;w<N;w++){
							field_on_target_support[Field_On_Support_Index(s,w,j)]
								=Scattered_Field_At_Point(X(m),Y(n),w,s,reflectivity,field_on_target_support,SP);
						}
					}
					j++;
				}
			}
		}
		std::cout<<"Generated "<<k<<" time Scattered field on target support"<<endl;

		//Generate data for each pulse
		for(int s=0;s<N;s++){
			//Generate data for each frequency
			for(int w=0;w<N;w++){
				//Find scattered field along flight path
				data[Data_Index(s,w)]+=Scattered_Field_At_Point(Gamma_X(S(s)),Gamma_Y(S(s)),w,s,reflectivity,field_on_target_support,SP);
				
				
			}
		}
		std::cout<<"Generated "<<k<<" time Scattered field in data"<<endl;
	}

	//Final step: Returns from the last scattering event to reach the antenna.
	//Generate data for each pulse
	for(int s=0;s<N;s++){
		std::cout<<"Working on slow time "<<s<<endl;
		//Generate data for each frequency
		for(int w=0;w<N;w++){
			//Find scattered field along flight path
			data[Data_Index(s,w)]+=Scattered_Field_At_Point(Gamma_X(s),Gamma_Y(s),w,s,reflectivity,field_on_target_support,SP);
		}
	}
	return data;
}


//Indexing
inline int Field_On_Support_Index(int s, int w, int j)
{

	return N*N*j+N*s+w;
}







