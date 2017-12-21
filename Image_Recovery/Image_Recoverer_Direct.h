//*
// Recover the target using direct integration of fourier data.
// This is a frequency domain method.
//
//
//
//
//
//*/
class ImageRecovererDirect : ImageRecoverer {};

void Image_Recovery_Direct(double *recovered_reflectivity,
                           complex<double> *data, ostream &log) {
  double A = -3200 * Band * pi / (M * c_0 * c_0);
  complex<double> *SUM = new complex<double>[M];
  for (int m = 0; m < N; m++) {
    for (int n = 0; n < N; n++) {
      SUM[Image_Index(m, n)] = 0;
    }
  }
  log << "Initialized sum" << endl;
  for (int m = 0; m < N; m++) {
    log << "Working on Row " << m << endl;
    for (int n = 0; n < N; n++) {

      for (int s = 0; s < N; s++) {
        for (int w = 0; w < N; w++) {
          SUM[Image_Index(m, n)] +=
              exp(-i * W(w) *
                  (2 * (X(m) * Gamma_D_X(S(s)) + Y(n) * Gamma_D_Y(S(s)) - r)) /
                  c_0) *
              conj(p(W(w))) * data[Data_Index(s, w)] /
              (W(w) * abs(p(W(w))) * abs(p(W(w))));
        }
      }
    }
  }

  log << "Sum  was computed" << endl;
  for (int m = 0; m < N; m++) {
    for (int n = 0; n < N; n++) {
      log << "Computed sum is:" << SUM[Image_Index(m, n)] << endl;
      recovered_reflectivity[Image_Index(m, n)] =
          A * abs(SUM[Image_Index(m, n)]);
    }
  }
  log << "Reflectivity was built" << endl;
}
}
;
