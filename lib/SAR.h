//*Represents a physical synthetic aperture radar system
//A SAR system has a time varying position. The system
//Emmits a time varying signal.  The system also measures
//the value of EMF on its antenna by some averaging process.
//To simplify things, we will model using the scalar wave equation. 
//All measurements are assumed to 
//
//
//This class should minimize the separation from software and hardware. 
//Ideally, we should be able to fill in an instance of SAR using a physical
//system.
//This system is a pulsed.  That is, it repeatedly sends out the signal p(t) 
//with frequency equal to the pulse rate. 
//*/

class SAR{
    public:
        Antenna antenna;

    Signal p;
    Carrier w0;
    float pulse_rate;
    float sample_rate;
    
    Data data;
    
    //The region in space where we will detect targets.  
    Scene scene;


    virtual point flight_path(double t);
    
};

///Represents the signal that the antenna will modulate and apply as a current to the 
//antenna. Assumes that the signal starts at t=ts, and ends at t=tf  
class Signal{

    float ts;
    float tf;

    virtual float p(float t);
    // Assume the samples are uniform and given at $(i+1/2N,j+1/2N)$ for
  // $i,j=0,1,\ldots N-1$  N is the number of slow time samples, frequency
  // samples. There are $N^2$ space samples
  const int N = 4;
 


};

class Chirp: Signal{
    float chirp_length;
    float start_frequency;
    
    
    complex<double> p(double w);{
	if((w<w_0+Band/2)&&(w>w_0-Band/2))
		return 1;
	else 
		return 0;
    }
    
};




class Rect: Signal{
    
    float p(float t){
        if(ts<t and t<tf )
            return 1;
        else return 0;
    }

};


///Antenna can act as tranmitter or reciever.
///As tranmitter, it takes the signal and creates an incident field.
///As reciever, it integrates the field incident on the antenna support.
///
class Antenna{
    ///The antenna center relative to the position of object it is mounted on. 
    ///center would be (0,0,0) if the antenna is at the center of the object. 
    point2d center;
    ///The antenna gain in a given direction
    double gain(vector3 direction);
    ///The 
    double recieve(signal p);
    double transmit(signal p);
};









}
	

