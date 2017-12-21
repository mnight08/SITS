#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

/// Take a given sar system and convert its collected data into an estimate for
/// the object that create that data.  This is an abstract class. Actually
/// recoverying and image depends on the geometry of the scenario and collected
/// data.
class ImageRecoverer {
  Target estimate;
  /// The system that collected the data
  SAR sar;

  virtual void write(string filename, string format);
  virtual void recover_target;

  // number of scattering events to consider in image recovery.
  const int K = 1;

  // image plane height
  const double image_height = 0;
};

/// The plane that
class ImagingPlane {
  const int M = N * N;
};

class ImageRecoverer2d : ImageRecoverer {};

class ImageRecoverer3d : ImageRecoverer {};
