#include "input.h"
//for 1d grid
inline int Grid_Index(int Box,int t, int l,ostream &log){
	return ((int) pow(2.0,(double) l)-1)*q+Box*q+t;
}


inline point Grid_Point(double * grids, int Box, int t1, int t2, int l,ostream &log){
	int B1=Box/N;
	int B2=Box%N;
    return point(grids[Grid_Index(B1,t1,l)] ,grids[Grid_Index(B2,t2,l)]);

}

//give the chebyshev grid for [0,1] with q points, sorted from left to right.  
inline void Proto_Grid(double *fillme,ostream &log){
	for(int t=0;t<q;t++){
		fillme[t]=.5*cos((q-1-t)*pi/(q-1))+.5;
	}
}

//fill in the array with the chebyshev grid adapted to box  at level l.
//Assumes that protogrid is already allocated.  b runs from 0 to 2^level-1
inline void Grid(double *grid, double *proto_grid, int Box, int l,ostream &log){
	double scale=pow(2,(double)-l);
	for(int t=0;t<q;t++){
		//shift protogrid by scale, to  grid from 0 to b.  Then shift by b*scale for the 

		//corresponding box.
		grid[t]=scale*(proto_grid[t]+Box);
	}
}

//Fills grids with all the 1d grids for [0,1] we will need. There are 2N-1 grids, and each has q entries.
//The function Grid_Index returns the appropriate index.  A tensor grid is used for higher dimensions
inline void Build_Grids(double* &grids,ostream &log){
	grids= new double[(2*N-1)*q];
	//make the first q entries the chebyshev grid adapted to [0,1].
	double * current_grid=grids+q;
	Proto_Grid(grids,log);
	//go through each scale
	for( int l=1;l<=L;l++){
		//go through each offset at that scale.
		for(int Box=0;Box<(int)(pow(2.0,l));Box++){
			Grid(current_grid,grids,Box,l);
			current_grid=current_grid+q;
		}
	}

}


inline void Log_Grids(double* &grids,ostream &log){

}

//The t-th Lagrange basis polynomial for the grid evaluated at the point x
inline double L1d(double *grids, int l, int Box, int t,double x,ostream &log){
	double product=1;
	int grid=Grid_Index(Box,q,l);
	log<<"T!!"<<t;
	for(int j=0;(j<q)&&(j!=t);j++){
		product*=(x-grids[grid+t])/(grids[grid+t]-grids[grid+j]);
	}
	return product;
}

//Tensor interpolation L_t1(y.x)*L_t2(y.y). 
//t goes up to q^2. B goes to N^2, so that t1, t2 are between 0, and q-1,
//and B1, and B2 are between 0, and M-1
inline double L2d(double *grids, int t1, int t2, int Box, int l,point y,ostream &log){
    int B1=Box/N;
	int B2=Box%M;
    return L1d(grids,l,B1,t1,y.x)*L1d(grids,l,B2,t2,y.y);
}



inline void Read_In_Data(complex<double> *data,string file,ostream &log){
	ifstream input;
	input.open(file.c_str());
	for(int s=0;s<N;s++) {
		for(int w=0;w<N;w++){
		   input>>data[Data_Index(s,w)];
		}
	}
	input.close();
}

inline void Swap(complex<double>* &current, complex<double>* &previous,ostream &log){
	complex<double>* temp;
	temp=previous;
	previous=current;
	current=temp;	
}

inline void Zero(complex<double> * weights ,ostream &log){
	for(int n=0;n<M*r_eps;n++)weights[n]=0;
}

//return the center of the box from the left at level l. 
inline double Center1d(int l, int box,ostream &log){
	double scale=pow(2,(double)-l);
	return (int)(scale*((double)box+0.5));
}

inline point Center(int l, int box,ostream &log){
	int b1=box/N;
	int b2=box%N;
	return point(Center1d(l,b1),Center1d(l,b2));
}

//We need to enumerate the boxes, and pairings we are working with.  
//Consistency is key.  In X_l there are 4^l cells.  In Y_l there are
//4^{L-l}.  The number of pairs is 4^L=N^2 for any level. 
//Each pair (A,B) at a scale corresponds to a number n=0,...N^2-1.
//A=0,...,4^l-1, and B=0,...,4^{L-l}-1.  A simple map
//is to take the first 2l bits of n to be A, and the last 2(L-l)
//bits to be B n->(n/4^l-1,4n^l)=(A,B).
//We always count from left to right, bottom up, so that we would start around (0,0), and end around
//(1,1) 
inline int A(int n, int l,ostream &log){
	return (n/(int)pow(2.0,2*(L-l)));
}

//shift n to the left, so that only the last 2*(L-l) digits are non Zero, 
//then shift back so that the first 2l bits are Zero.
inline int B(int n, int l,ostream &log){
	return n%(int)(pow(2.0,2*(L-l)));
}

//Parent of A is just A with its last two bits dropped
inline int Parent(int A,ostream &log){
	return A>>2;
}

//Children of B are numbers where we append two bits
inline int Child(int B,int c,ostream &log){
    return (B<<2)+c;
}
//There are N*N*r_eps weights to index.  Here, B needs to have at most 2l bits 
inline int Weight_Index(int A, int B, int l, int t1, int t2,ostream &log){
	return ((A<<(2*(L-l)))+B)*r_eps+q*t1+t2;
}

///Find the box that p is in  at level l. l=0 means only one cell.  
//The unit box is naturally partitioned. The box that a point is in is completely determined by
//the first 2l bits of its components. First l (from left) bits are where y is, next l are x.
inline int Box(point p, int l,ostream &log){
	//shifts x, and y left by l bits by multiplying by 2^l.
	int column =(int) (p.x*pow(2.0,(double) l));
	int row    =(int) (p.y*pow(2.0,(double) l));
	return row*pow(2.0,(double) l)+column;
}
inline int Image_Index(int m, int n){
	return N*n+m;
}
inline int Data_Index(int s,int w){
	return N*s+w;
}
inline double X(int m){
	return m/(double)N+ 1.0/(2.0*N);
}
inline double Y(int n){
	return n/(double)N+ 1.0/(2.0*N);
}
inline double Gamma_D_X(double sl){
	return Gamma_X(sl)/sqrt(Gamma_X(sl)*Gamma_X(sl)+Gamma_Y(sl)*Gamma_Y(sl));
}

inline double Gamma_D_Y(double sl){
	return Gamma_Y(sl)/sqrt(Gamma_X(sl)*Gamma_X(sl)
					+Gamma_Y(sl)*Gamma_Y(sl));
}
inline double W(int w){
	return w_0-Band/2+w*Band/(double)N+1.0/(2.0*N);
}

inline double Q(int e){
		  return e/N+1.0/(2*N);
}
inline double S(int s){
	return s/(double)N+1.0/(2*N);
}

inline void Zero_Grids(double *grids){
	for(int i=0;i<((int)pow(2.0,L+1)-1)*q;i++){
		grids[i]=0.0;
	}
}


//Perform image recovery using butterfly algorithm
//Parameters are defined in input.h. Conventions are handled by functions above.  
//We follow the algorithm developed in: "A butterfly algorithm for synthetic aperture radar imaging"
//At each level there are N^2q^2 weights to compute. We need weights from the previous level
//to compute weights at current level.
void Image_Recovery_Butterfly(double*& recovered_reflectivity, complex<double>* data,ostream &log){
	int l=0;
	//Constants outside the sum
	double C=-3200*Band*pi/(M*c_0*c_0),qj=0,sl=0;
	
	point x,y;
	//complex valued sum. We will need to take modulus before
	//it is reflectivity estimate.
	complex<double> *SUM=new complex<double>[M];
	recovered_reflectivity=new double[M];
	
	for(int n=0;n<N;n++){
		for(int m=0;m<N;m++){
			SUM[Image_Index(m,n)]=0;
			recovered_reflectivity[Image_Index(m,n)]=0;
		}
	}

	complex<double> *previous_weights=new complex<double>[M*r_eps];
	complex<double> *weights = new complex<double>[M*r_eps];

	log<<"Made it here"<<endl;
	//Compute all the 1-d grids we will need.
	double *grids;
	Build_Grids(grids);
	Print_Grids(grids);
	log<<"Well the grids were built"<<endl;
	//Step 1: Initialize the weights. Done different than the paper to not need tree.
	//go through all the points that we are integrating over and add them to the weight
	//for the box they are in. This is so we do not need to have an actual tree structure
	//The exponential is moved inside the sum.
  	Zero(weights);
	Zero(previous_weights);
 	for(int e=0;e<N;e++){
		for(int o=0;o<N;o++){
			y=point(S(o),Q(e));
			if(abs(f(y,data))!=0){
				for(int t1=0;t1<q;t1++){
					for(int t2=0;t2<q;t2++){
						weights[Weight_Index(0,Box(y,L),l,t1,t2)]+=
							exp(-i*(phi(Center(0,0),Grid_Point(grids,Box(y,L),t1,t2,l))
								-phi(Center(0,0),y)))
							*L2d(grids,t1,t2,Box(y,L),L,y)
							*f(y,data);
					}
				}
			}
		}
	}


	
	log<<"First Level is done"<<endl;
	Swap(previous_weights,weights);
	Zero(weights);
	//Step 2: for l=1,..., L/2
	//Compute the next weights from previous weights using the fact when diam(B)<1/N^2
	for(l=1;l<=L/2;l++){
		log<<"Current level is"<<l<<endl;
		//go through all the pairs of boxes.  There are N^2.
		for(int n=0;n<M;n++){
			//go through the children of B paired with the parent of A to compute weights
			for(int t1=0;t1<q;t1++){
				for(int t2=0;t2<q;t2++){
					for(int c=0;c<4;c++){
						//for each child parent pair
						for(int tc1=0;tc1<q;tc1++){
							for(int tc2=0;tc2<q;tc2++){
								if(abs(previous_weights[Weight_Index(Parent(A(n,l)),Child(B(n,l),c),l-1,tc1,tc2)])!=0){
									weights[Weight_Index(A(n,l),B(n,l),l,t1,t2)]+=
											  L2d(grids, tc1,tc2, B(n,l),L-l,Grid_Point(grids,Child(B(n,l),c),tc1,tc2,L-sl+1))
									*exp(i*phi(Center(A(n,l),l),Grid_Point(grids,Child(B(n,l),c),tc1,tc2,L-l+1)))
									*previous_weights[Weight_Index(Parent(A(n,l)),Child(B(n,l),c),l-1,tc1,tc2)];
								}	
							}
						}	
					}
					weights[Weight_Index(A(n,l),B(n,l),l,t1,t2)]
							  *=exp(-i*phi(Center(A(n,l),l),Grid_Point(grids,B(n,l),t1,t2,L-l)));
				}
			}
		}
		Swap(previous_weights,weights);
		Zero(weights);
	}


	log<<"Step 2 done!!"<<endl;
	//Step 3: for l=L/2 Do switch
	for(int n=0;n<M;n++){
		for(int t1=0;t1<q;t1++){
			for(int t2=0;t2<q;t2++){
				for(int s1=0;s1<q;s1++){
					for(int s2=0;s2<q;s2++){
						if(abs(previous_weights[Weight_Index(A(n,L/2),B(n,L/2),L/2,s1,s2)])!=0){
							weights[Weight_Index(A(n,L/2),B(n,L/2),L/2,t1,t2)]
								+=exp(i*phi(Grid_Point(grids,A(n,l),t1,t2,l),Grid_Point(grids,B(n,l),s1,s2,l)))
								*previous_weights[Weight_Index(A(n,L/2),B(n,L/2),L/2,t1,t2)];
						}
					}

				}
			}
		}
	}
	log<<"Step 3 is done!!"<<endl;
	Swap(previous_weights,weights);
	Zero(weights);
	//Step 4: for l=L/2,...,L Compute the next weights from previous 
	//weights using the fact that diam(A)<1/M handle l=L/2+1,... L
  	for(;l<=L;l++){
		log<<"Current level is"<<l<<endl;
		for(int n=0;n<M;n++){
			//set grid for current box b.
			for(int t1=0;t1<q;t1++){
				for(int t2=0;t2<q;t2++){
					for(int c=0;c<4;c++){
						//for each child parent pair
						for(int tp1=0;tp1<q;tp1++){
							for(int tp2=0;tp2<q;tp2++){
								//The indices for the pairings are given by
								if(abs(previous_weights[Weight_Index(Parent(A(n,l)),Child(B(n,l),c),l-1,tp1,tp2)])!=0){
									weights[Weight_Index(A(n,l),B(n,l),l,t1,t2)]+=L2d(grids,tp1,tp2,Parent(A(n,l)),l-1,Grid_Point(grids,A(n,l),t1,t2,l))
									*exp(-i*phi(Grid_Point(grids,Parent(A(n,l)),tp1,tp2,L-l+1),Center(L-l+1,Child(B(n,l),c))))
									*previous_weights[Weight_Index(Parent(A(n,l)),Child(B(n,l),c),l-1,tp1,tp2)];
								}
							}
						}
						weights[Weight_Index(A(n,l),B(n,l),l,t1,t2)]
							=exp(i*phi(Grid_Point(grids,A(n,l),t1,t2,l),Center(l,Child(B(n,l),c))))
							*previous_weights[Weight_Index(Parent(A(n,l)),Child(B(n,l),c),l,t1,t2)];
					}
				}
			}
		}
		Swap(previous_weights,weights);
		Zero(weights);
	}
	log<<"Step 4 is done!!!!!!!!!!!!!!!!"<<endl;
	//undo the last Swap since we need the latest computed weights.
	Swap(previous_weights,weights);

	//Step 5: use equivalent weights to compute full sum. B=0 since we are doing the entire integral
	for(int n=0;n<N;n++){
		for(int m=0;m<N;m++){
			x=point(X(m),Y(n));
			//Get the box A that x is in
			for(int t1=0;t1<q;t1++){
				for(int t2=0;t2<q;t2++){
					SUM[Image_Index(m,n)]+=L2d(grids,t1,t2,Box(x,L),L,x)
							  *exp(-i*phi(Grid_Point(grids,0,t1,t2,0),Center(0,0)))
							  *weights[Weight_Index(Box(x,L),0,L,t1,t2)];
				}
			}
			SUM[Image_Index(m,n)]=C*exp(i*phi(x,Center(0,0)))*SUM[Image_Index(m,n)];
			recovered_reflectivity[Image_Index(m,n)]=abs(SUM[Image_Index(m,n)]);
		}
	}
	//clean up memory
	delete[] weights;
  	delete[] previous_weights;
	delete[] SUM;
}


void Image_Recovery_Direct(double * recovered_reflectivity, complex<double>* data, ostream &log){
	double A=-3200*Band*pi/(M*c_0*c_0);
	complex<double> * SUM= new complex<double>[M];
	for(int m=0;m<N;m++){
		for(int n=0;n<N;n++){
			SUM[Image_Index(m,n)]=0;	
		}
	}
	log<<"Initialized sum"<<endl;
	for(int m=0;m<N;m++){
		log<<"Working on Row "<<m<<endl;
		for(int n=0;n<N;n++){

			for(int s=0;s<N;s++){
				for(int w=0;w<N;w++){
					SUM[Image_Index(m,n)]+=exp(-i*W(w)*(2*(X(m)*Gamma_D_X(S(s))+Y(n)*Gamma_D_Y(S(s))-r))/c_0)
							  *conj(p(W(w)))*data[Data_Index(s,w)]/(W(w)*abs(p(W(w)))*abs(p(W(w))));	
				}
			}
		}	
	}

	log<<"Sum  was computed"<<endl;
	for(int m=0;m<N;m++){
		for(int n=0;n<N;n++){
			log<<"Computed sum is:"<<SUM[Image_Index(m,n)]<<endl;
			recovered_reflectivity[Image_Index(m,n)]=A*abs(SUM[Image_Index(m,n)]);	
		}
	}
	log<<"Reflectivity was built"<<endl;
}

//Output reconstruction for viewing
void Output_Reflectivity(double * reflectivity,string file){
	ofstream output;
	output.open(file.c_str());
	for(int n=0;n<N;n++){
		for(int m=0;m<N;m++){
		   output<<reflectivity[Image_Index(m,n)]<<" ";
		}
		output<<endl;
	}
	output.close();
}
