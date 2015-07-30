

#ifndef _LIB_IO_GAUSSCONV_H_
#define _LIB_IO_GAUSSCONV_H_

double* read_image_to_gray(char* fname, int *pw, int *ph);


void printImage(double* imageIn, int w, int h, char* name);                     

void printImage_LinearConversion(double *imageIn, int w, int h, char *name);



#endif /* _LIB_IO_GAUSSCONV_H_ */




