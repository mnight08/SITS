const complex<double> i(0, 1);
const double pi = 3.14159265359;
const double c_0 = 299792458;

class Math {

  inline double Norm(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
  }

  // Linear interpolation.  Assume that data is sorted in frequency and slow
  // time.
  complex<double> interp(complex<double> *d, double sl, double qj) {
    int s = (int)(sl * N);
    int w = (int)(qj * N);
    return d[Data_Index(s, w)];
  }
class point{
	public:
	double x, y;
	inline void set(double xp,double yp){x=xp;y=yp;}
	point(){x=0;y=0;}
	point(double xp, double yp){x=xp;y=yp;}
};

// The t-th Lagrange basis polynomial for the grid evaluated at the point x
  inline double L1d(Grid grid, int l, int Box, int t, double x, ostream &log) {
    double product = 1;
    ChebyshevGrid grid = ChebyshevGrid(Box, q, l);
    log << "T!!" << t;
    for (int j = 0; (j < q) && (j != t); j++) {
      product *=
          (x - grid.values[grid + t]) / (grids[grid + t] - grids[grid + j]);
    }
    return product;
  }

  // Tensor interpolation L_t1(y.x)*L_t2(y.y).
  // t goes up to q^2. B goes to N^2, so that t1, t2 are between 0, and q-1,
  // and B1, and B2 are between 0, and M-1
  inline double L2d(double *grids, int t1, int t2, int Box, int l, point y,
                    ostream &log) {
    int B1 = Box / N;
    int B2 = Box % M;
    return L1d(grids, l, B1, t1, y.x) * L1d(grids, l, B2, t2, y.y);
  }
};
}

