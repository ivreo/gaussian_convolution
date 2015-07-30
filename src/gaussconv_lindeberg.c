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
    fprintf(stderr, "Lindeberg convolution  \n");
    fprintf(stderr, "     gaussconv_lindeberg in sigma gamma(laplacian scheme) out  \n");
}



int main(int argc, char **argv)
{
    double sigma;
    double g;
    switch (argc){
        case 5:  
            sigma = atof(argv[2]);
            g = atof(argv[3]);
            break;
        default:
            print_usage();
            return -1;
    }
    int w, h;
    double* x = read_image_to_gray(argv[1], &w, &h);
    double* y = (double*)malloc(w*h*sizeof(double));
    gaussconv_lindeberg(x, y, w, h, sigma, g);

    printImage(y, w, h, argv[4]);
    free(x);
    free(y);
    return 0;
}
