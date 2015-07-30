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
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lib_io_gaussconv.h"
#include "lib_gaussconv.h"





void print_usage()
{
    fprintf(stderr, "DEMO semigroup Lindeberg convolution  \n");
    fprintf(stderr, " gaussconv_lindeberg in sigma gamma(laplacian scheme) N(iterations) out  \n");
}


double rmse(double* x1, double* x2, int l)
{
    double diff;
    double rmse = 0.0;
    for(int i = 0; i < l; i++){
        diff = x1[i]-x2[i];
        rmse += diff*diff;
    }
    rmse = sqrt( rmse / (double)l )/255.;
    return rmse;
}




int main(int argc, char **argv)
{
    double sigma;
    double g;
    int N;

    switch (argc){
        case 6:  
            sigma = atof(argv[2]);
            g = atof(argv[3]);
            N = atoi(argv[4]);
            break;
        default:
            print_usage();
            return -1;
    }
    int w, h;
    double* x = read_image_to_gray(argv[1], &w, &h);
    double* t = (double*)malloc(w*h*sizeof(double));  //temporary
    double* y1 = (double*)malloc(w*h*sizeof(double)); //direct
    double* y2 = (double*)malloc(w*h*sizeof(double)); //iterated
    double* d = (double*)malloc(w*h*sizeof(double));  //difference

    // direct
    gaussconv_lindeberg(x, y1, w, h, sqrt((double)N)*sigma, g);

    // iterated
    for(int j=0; j<w*h; j++)
        t[j] = x[j];
    for(int i=0; i<N;i++){
        gaussconv_lindeberg(t, y2, w, h, sigma, g);
        for(int j=0; j<w*h;j++)
            t[j] = y2[j];
    }

    // difference image
    for(int j=0; j<w*h; j++){
        d[j] = y1[j]-y2[j];
    }

    // output 
    char name[FILENAME_MAX];
    
    sprintf(name, "%s_direct.png", argv[5]);
    printImage(y1, w, h, name);

    sprintf(name, "%s_iter.png", argv[5]);
    printImage(y2, w, h, name);
    
    sprintf(name, "%s_diff.png", argv[5]);
    printImage_LinearConversion(d, w, h, name);

    double err = rmse(y1, y2, w*h);
    fprintf(stdout, "RMSE = %g \n", err);
    
    free(x);
    free(t);
    free(y1);
    free(y2);
    free(d);

   
    return 0;
}
