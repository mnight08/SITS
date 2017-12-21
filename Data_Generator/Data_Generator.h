const complex<double> i(0,1);
const double pi=3.14159265359;
const double c_0=299792458;
//Assume the samples are uniform and given at $(i+1/2N,j+1/2N)$ for $i,j=0,1,\ldots N-1$
//N is the number of slow time, frequency samples.  There are $N^2$ space samples
const int N=4;
const int M=N*N;
//number of scattering events to consider.


const int K=1;
//image plane height
const double image_height=0;


//normalized function domain
complex<double> f(point y,complex<double>*d){
	return (conj(p(w_0-Band/2+y.y*Band))*interp(d,y.x,y.y))
		/((w_0-Band/2+y.y*Band)*abs(p(w_0-Band/2+y.y*Band))*abs(p(w_0-Band/2+y.y*Band)));
}

class DataGenerator{
    //N is number of bounces.
    virtual void GenerateData(float [] x, float [] y, n, m, int N);
    virtual void LoadTarget(Target& T);
    virtual void Write();

    Target target;  
    Data data;

    DataGenerator{
    
    }

};

class DataGenerator2d: DataGenerator{


    void write(ofstream &file)
    {

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
    }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~///




}


///The data that we would collect physically.  Data may have several dimensions.
///For example, there are three polarizations, slow time, fast time, flight pass, etc.
///Data can eithe be representedin frequency domain, or time domain.  Typically if 
///data is given in frequency domain, then for each pulse(slowtime bin), The collected 
///data has a windowed Fourier transform aplied to it. This is either done in hardward, 
///digitally.
class Data{
    String format="RAW";

};

///
class Data2d:Data{

  void load_data(string file) {
    ifstream input;
    input.open(file.c_str());
    for (int s = 0; s < N; s++) {
      for (int w = 0; w < N; w++) {
        input >> data[Data_Index(s, w)];
      }
    }
    input.close();
  

}

/**
 *A SAR is a system that has a moving antenna which coherently 
 *
 *
 * **/
class SAR{
    point2d flight_path(float t);
    point2d 
//center frequency
const double w_0=pow(10,9);
//Bandwidth
const double Band=pow(10,8);
//radius of flight path
const double r=100;
//Flight height
const double H=100;

}

class DataManifold{
    virtual write();
};

class Field{

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


///Generate data due to the target with number of scattering events K
void Generate_Data2dNaive(int K)
{
    //Data is complex valued.  There are N^2=M data points.
    //complex<double> *data=new complex<double>[M];
    complex<double> *data=new complex<double>[N*N];
  		
    int SP=target.Size_Of_Support();
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

};




 
void Generate_Data(string infilename){
    target.load(infilename);
    data=Generate_Data(reflectivity,1);
    string outfilename=infilename.erase(infilename.end()-3,infilename.end())+".out";
    Output_Data(data,outfilename);
}








