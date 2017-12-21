#include "image_recovery.h"
#include "input.h"

// This file will test the ability to recover reflectivities. The input is
// complex valued data file.  The output is an esimate to the real valued
// reflectivity bitmap
int main(int argc = 1, char *filenames[]) {

  Data data(filenames);

  // The reflectivity will be a list of $N^2$ real numbers.
  Target recovered_target;
  Read_In_Data(data, "../Data/" + filenames + ".in");

  Image_Recovery_Direct(recovered_reflectivity, data);
  Output_Reflectivity(recovered_reflectivity,
                      "../Recovered_Targets/" + filenames + "_Direct.out");

  Image_Recovery_Butterfly(recovered_reflectivity, data);
  Output_Reflectivity(recovered_reflectivity,
                      "../Recovered_Targets/" + filenames + "_Butterfly.out");
  return 0;
}
