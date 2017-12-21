/**A target represents a physical entity that scatters electro magnetic radiation.
 *In the wave model, any target is embedded in the reflectivity function. The reflectivty function
 *is represented as a space and time dependent real valued function whose value is proportional to 
 *the shinyness of the target. Currently, targets are assumed stationary. A target should be able to
 *save itself to file in approriate format, load from file, 
 *Subclasses include sampled targets, and function
 *targets.  A function target takes in a point and returns the reflectivity at that point.  A sampled target
 *is a discrete representation of the physical entity we are modeling.  
 *
 *
 * */


class Target{
    virtual operator () (double x, double y);
    virtual void write();
    virtual void load();
    virtual double reflectivity();

    String format;
        
};


///A sampled 2d target, consist of a collection of samples of the reflectivity function.
///These samples are used to evaluate teh relfectivity function using some interpolation method.
class Target2dS:Target{
    double [][] samples;

    //Number of non zero entries in $v$
int Size_Of_Support(double *reflectivity){
	int size=0;
	for(int m=0;m<N;m++)
		for(int n=0;n<N;n++)
			if(reflectivity[Image_Index(m,n)]!=0)
				size++;
	return size;
}
    void Load(string name)
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

    };

class Target3dS:Target{};

class Target2dF:Target{
    
};



