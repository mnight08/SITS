/**
 *A target represents a physical entity that we would scatter
 an incident field. Curently, targets are assumed stationary. Should relax this latter.  Subclasses include sampled targets, and function
 targets.  A function target takes in a point and returns the reflectivity at that point.  
 *
 *
 * */
class Target{
    virtual operator () (double x, double y);
    virtual void write();
    virtual void load()
    {
        
        
    }
    String format;
    double [][] reflectivity;

        
};
class point{
	public:
	double x, y;
	inline void set(double xp,double yp){x=xp;y=yp;}
	point(){x=0;y=0;}
	point(double xp, double yp){x=xp;y=yp;}
};
complex<double> f(point y,complex<double>*d);
complex<double> p(double w);
complex<double> interp(complex<double>* d, double sl,double qj);
double phi(point x,point y);
double Gamma_X(double s);
double Gamma_Y(double s);

class Target2dS{};

class Target3dS{};

class Target2dF{};


//Number of non zero entries in $v$
int Size_Of_Support(double *reflectivity){
	int size=0;
	for(int m=0;m<N;m++)
		for(int n=0;n<N;n++)
			if(reflectivity[Image_Index(m,n)]!=0)
				size++;
	return size;
}

    void LoadTarget(string name)
    {
    
        //Reflectivities are real.  There are N^2=M image points.  
            
        ifstream input;
        input.open((infile_name+".in").c_str());
        for(int n=0;n<N;n++) {
	    for(int m=0;m<N;m++){
	        input>>reflectivity[Image_Index(m,n)];
	    }
        }
        input.close();
        }       

    }
