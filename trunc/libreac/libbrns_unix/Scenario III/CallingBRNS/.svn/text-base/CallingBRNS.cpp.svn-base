// CallingBRNS.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>


#ifdef UNIX_GNU
#include <dlfcn.h>
#endif

#ifndef UNIX_GNU
extern "C" __declspec( dllimport )  void brnsIsAlive();
extern "C" __declspec( dllimport )  void invokebrns(double* theArray, double* theArray2, double* outputArray, int sizeOfArray, double timeStep, int* fixedConcBndry, int* returnValue, double pos_x, double pos_y, double pos_z);
#endif

using namespace std;

void printArray(double* arrayToPrint, int sizeOfArray)
{
     for (int i = 0 ; i < sizeOfArray; i++){
          cout << i << ": " << arrayToPrint[i] << endl;
     }
}

#ifndef UNIX_GNU
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char* argv[])
#endif
{
     #ifdef UNIX_GNU
     void *hDll = dlopen("./brns.so", RTLD_NOW);
     typedef void (* LPFNDLLFUNC)(double*, double*, double*, int&, double&, int*, int*);
     LPFNDLLFUNC invokebrns = NULL;
     invokebrns = (LPFNDLLFUNC)dlsym(hDll, "invokebrns_");
     #endif

     cout << "Hello World!" << endl;
     cout << "Trying to call BrnsDLL ..." << endl;
//     brnsIsAlive();
//     cout << "Hm, looks good, I guess." << endl;


     double* myArray = NULL;
     double* myArray2 = NULL;
     double* outputArray = NULL;
     int* fixedConcBndry = NULL;
     double timeStep;
     int sizeOfArray;
     int returnValue = -111;

     timeStep = 10000;

     sizeOfArray=6;
     myArray = new double[sizeOfArray];
     myArray2 = new double[sizeOfArray];
     outputArray = new double[sizeOfArray];
     fixedConcBndry = new int[sizeOfArray];

     myArray[0] = 0.12;
     myArray[1] = 0.34;
     myArray[2] = 0.56;
     myArray[3] = 0.78;
     myArray[4] = 0.910;
     myArray[5] = 0.919;
     myArray2[0] = 0.11;
     myArray2[1] = 0.33;
     myArray2[2] = 0.55;
     myArray2[3] = 0.77;
     myArray2[4] = 0.99;
     myArray2[5] = 0.991;
     outputArray[0] = -1.0;
     outputArray[1] = -1.0;
     outputArray[2] = -1.0;
     outputArray[3] = -1.0;
     outputArray[4] = -1.0;
     outputArray[5] = -1.0;
     fixedConcBndry[0] = 0;
     fixedConcBndry[1] = 0;
     fixedConcBndry[2] = 0;
     fixedConcBndry[3] = 0;
     fixedConcBndry[4] = 0;
     fixedConcBndry[5] = 0;


     cout << "Concentrations after Transport:" << endl;
     printArray(myArray, sizeOfArray);

     cout << "Calling invokebrns() ..." << endl;

     invokebrns(myArray, myArray2, outputArray, sizeOfArray, timeStep, fixedConcBndry, &returnValue, 1.1,2.2,3.3);

     cout << "Returning from invokebrns(), return code " << returnValue << endl;
     cout << "New Concentrationvector:" << endl;

     printArray(outputArray, sizeOfArray);
     
     delete [] myArray;
     delete [] myArray2;
     delete [] outputArray;
     delete [] fixedConcBndry;
     myArray = NULL;
     myArray2 = NULL;
     outputArray = NULL;
     fixedConcBndry = NULL;

     #ifdef UNIX_GNU
     dlclose(hDll);
     #endif
     return 0;
}
