
// Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
// // <ives.rey-otero@cmla.ens-cachan.fr>
// //
// // Version 20140515 (May 15th, 2014)
// //
// //
// // This program is free software: you can use, modify and/or
// // redistribute it under the terms of the GNU General Public
// // License as published by the Free Software Foundation, either
// // version 3 of the License, or (at your option) any later
// // version. You should have received a copy of this license along
// // this program. If not, see <http://www.gnu.org/licenses/>.
//
//
//


#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "io_png.h"

#include "lib_io_gaussconv.h"
#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )


double* read_image_to_gray(char* fname, int *pw, int *ph)
{
    size_t w, h;
    float* imf32 = io_png_read_f32_rgb(fname, &w, &h);
    double* x = (double*)malloc(w*h*sizeof(double));
    for(int i = 0; i < w*h; i++)
        x[i] = (double)((0.299*imf32[i] + 0.587*imf32[w*h+i] + 0.1140*imf32[2*w*h+i])); //RGB2GRAY
    free(imf32);
    *pw = w;
    *ph = h;
    return x;
}




void printImage(double *imageIn, int w, int h, char *name)
{
    float *imf32 = (float *) malloc(w * h * sizeof(float));
    for (int i = 0; i < w * h; i++) {
        imf32[i] =  (float) imageIn[i];
    }
    io_png_write_f32(name, imf32, w, h, 1);
    free(imf32);
}





void capped_linear_conversion(double *in, double *out, int l, double mindiff)
{
    double a, b;
    double minV = DBL_MAX;
    double maxV = DBL_MIN;
    for (int i = 0; i < l; i++) {
        if (in[i] >= maxV){
            maxV = in[i];
        }
        if (in[i] <= minV){
            minV = in[i];
        }
    }
    double div = MAX( maxV - minV, mindiff );
    a = (250.0 - 0.0) / div ;
    b = -a * minV;
    for (int i = 0; i < l; i++)
        out[i] = a * in[i] + b;
}



void linear_conversion(double *imIn, double *imOut, int length){
    double a, b;
    double minVal = DBL_MAX;
    double maxVal = DBL_MIN;
    for (int i = 0; i < length; i++) {
        if (imIn[i] >= maxVal){
            maxVal = imIn[i];
        }
        if (imIn[i] <= minVal){
            minVal = imIn[i];
        }
    }
    a = (250.0 - 0.0) / (maxVal - minVal);
    b = -a * minVal;
    for (int i = 0; i < length; i++)
        imOut[i] = a * imIn[i] + b;
}

void printImage_LinearConversion(double *imageIn, int w, int h, char *name){
    double *imTemp = (double*)malloc(w*h*sizeof(double));
    linear_conversion(imageIn, imTemp, w * h);
    //capped_linear_conversion(imageIn, imTemp, w * h, 0.00000001);
    linear_conversion(imageIn, imTemp, w * h);
    printImage(imTemp, w, h, name);
    free(imTemp);
}


