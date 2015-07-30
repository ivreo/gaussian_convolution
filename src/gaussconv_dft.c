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

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lib_io_gaussconv.h"
#include "lib_gaussconv.h"



void print_usage()
{
    fprintf(stderr, "smoothing via DFT convolution  \n");
    fprintf(stderr, "   gaussconv_dft IN SIGMA OUT  \n");
}



int main(int argc, char **argv)
{
    double sigma;
    switch (argc){
        case 4:
            sigma = atof(argv[2]);
            break;
        default:
            print_usage();
            return -1;
    }
    int w, h;
    double* x = read_image_to_gray(argv[1], &w, &h);
    double* y = (double*)malloc(w*h*sizeof(double));
    gaussconv_dft(x, y, w, h, sigma);
    
    printImage(y, w, h, argv[3]);
    free(x);
    free(y);
    return 0;
}
