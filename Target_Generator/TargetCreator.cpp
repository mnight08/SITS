#include <fstream>
#include <stdio.h>      /* printf */
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;
//useage: TargetCreator [filename]
//If no argumetns, then just open paint, if there is a file,
//then open that file with paint.  In the future paint should be replaced
//with some non proprietary image editor.


int main(int arg_count,char *arguments[]){
    //if the creator has no arguments, then open paint for now.  Will replace
    //if it is useful.
//
    //cout<<arg_count<<endl;
    //cout<<arguments;
    //first argument is always the path
    if(arg_count==1)  
        system("mspaint");
    //Second argument will be a file to open
    else if (arg_count==2)
    {
        string path="../Input/Targets/";
        //for( int i=0;i<2;i++)
        //{
        //    cout<<string(arguments[i])<<endl;
        //}
       
        cout<<"mspaint "+path+string(arguments[1]);
        system(("mspaint "+path+string(arguments[1])).c_str());
    }


   //fstream savefile.openfile("target_editor.conf");
   // fstream savefile.openfile("target_conf");




}

//void edit_target(fstream &file){
    

//}

//void start_editor();
