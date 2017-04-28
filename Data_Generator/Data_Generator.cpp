#include "Data_Generator.h"

//This file will generate data and store it in the data folder
//Matlab will be used for reading in a bitmap representing the reflectivity and converting 
//it to a list of numbers. This program will use the list to generate complex valued data 
//and store it in a file for use by the Image_Recovery program.
//Input is a filename to get reflectivity from and a file name to output the data to.
//this assumes that the reflectivity functions is defined on a plane.
int main(int argc, char *argv[])
{
    if(arc>=2){
        for(int i=1;i<argc;i++){
            Generate_Data(argv[i]);
        }
    }
    return 0;
}
