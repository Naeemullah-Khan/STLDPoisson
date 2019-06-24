#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include "CG.h"

__inline int mindex3(int x, int y, int z, int sizx, int sizy) { return x+(y*sizx)+(z*sizx*sizy); }

template<class ImageType>
int WriteVector(std::string file_name, ImageType *vector_data, int length)
{
	FILE* file = fopen(file_name.c_str(), "wb");
	if (file == NULL)
	{
		std::cout << "Error: file for writing vector data could not be open!" << std::endl;
		return -1;
	}
	int image_size = length;
	fwrite(vector_data, sizeof(ImageType), image_size, file);
	fclose(file);
	return 1;
}

template <typename T>
std::string to_string2(T value)
{
	std::ostringstream os ;
	os << value ;
	return os.str() ;
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ) 
{
    
    // Check for proper number of input and output arguments
    if(nrhs<3) mexErrMsgTxt("3 inputs are required.");


    // Check data input types
    if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) mexErrMsgTxt("Images vector must be of class double");
    if(mxGetClassID(prhs[1])!=mxINT32_CLASS) mexErrMsgTxt("Alpha must be of class int");
    if(mxGetClassID(prhs[2])!=mxINT32_CLASS)  mexErrMsgTxt("Mask must be of class int");


    // Get the sizes of the input image volume
    const mwSize *dims;
    int height=0, width=0, feature_vector_size=0;
    if(mxGetNumberOfDimensions(prhs[0])==3) 
    {
        dims= mxGetDimensions(prhs[0]);
        height=dims[0];
        width=dims[1];
        feature_vector_size=dims[2];
    }
    else 
    {
        mexErrMsgTxt("Image must be 3d.");        
    }

    // read input
    double *I_big_vector_all=(double*)mxGetPr(prhs[0]);
    int *alpha_vector=(int*)mxGetPr(prhs[1]);
    int* mask=(int*)mxGetPr(prhs[2]);
	double *tol = (double*)mxGetPr(prhs[3]);
      
    // output
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *u_all= mxGetPr(plhs[0]);   
    int IJK_index;
    double *I;
    char* mask_in=new char[height*width];
    for(int i=0;i<height*width;i++)
    {
        mask_in[i]=(char)mask[i];
    }   
    //WriteVector("mask_in.dat", mask_in, height*width);
    //WriteVector("edgebits_in.dat", edgebits_in, height*width);
    //WriteVector("alpha_vector.dat", alpha_vector, feature_vector_size);
    for(int k=0;k<feature_vector_size;k++)
    {
        IJK_index=mindex3(0, 0, k, dims[0], dims[1]);
        //WriteVector("I_"+to_string2(k)+".dat", I_big_vector_all+IJK_index, height*width);
        ComputeCGSolution(tol[0], u_all+IJK_index, I_big_vector_all+IJK_index, mask_in, alpha_vector[k], height, width);
        //WriteVector("u_"+to_string2(k)+".dat", u_all+IJK_index, height*width);
        //printf("Thread: %d/%d\n", omp_get_thread_num(), omp_get_num_threads());
    }
 
    delete[] mask_in;
}