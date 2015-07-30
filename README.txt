Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20150415 (May 15th, 2014)


===============================================================================
== Overview ===================================================================

This C ANSI source code is linked to the following publication:

    [1] "Computing an exact Gaussian scale-space"
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2014.
        http://www.ipol.im/pub/algo/rd_semigroup/

A demonstration facility is available at:
        http://www.ipol.im/pub/demo/rd_semigroup/



===============================================================================
== Licence =================================================

This program is free software: you can use, modify and/or
redistribute it under the terms of the GNU General Public
License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later
version. You should have received a copy of this license along
this program. If not, see <http://www.gnu.org/licenses/>.



===============================================================================
== How to use this code =======================================================

-------------
----- 1) To apply a Gaussian smoothing of parameter SIGMA to an image
-------------

- with DCT convolution:
    gaussconv_dct IMAGE_IN SIGMA IMAGE_OUT

- with DFT convolution:
    gaussconv_dft IMAGE_IN SIGMA IMAGE_OUT

- with Sampled Gaussian kernel convolution
    gaussconv_sampled_kernel IMAGE_IN SIGMA k(truncation) IMAGE_OUT

- with Lindeberg convolution:
    gaussconv_Lindeberg IMAGE_IN SIGMA gamma(Laplacian scheme) IMAGE_OUT

For each method, the output image is IMAGE_OUT (png).

------------
----- 2) To check the semi-group property
         (i.e., To compare the result of one smoothing of
         parameter \sqrt(N) SIGMA with the result of N consecutive smoothings
         of parameter SIGMA)
-------------

- with DCT convolution:
    demo_gaussconv_dct IMAGE_IN SIGMA  N  LABEL_OUT

- with DFT convolution:
    demo_gaussconv_dft IMAGE_IN SIGMA N  LABEL_OUT

- with Sampled Gaussian kernel convolution
    demo_gaussconv_sampled_kernel IMAGE_IN SIGMA k(truncation)  N  LABEL_OUT

- with Lindeberg convolution:
    demo_gaussconv_Lindeberg IMAGE_IN SIGMA gamma(Laplacian scheme) N LABEL_OUT

For each method, the output images are
      LABEL_OUT_direct.png
      LABEL_OUT_semi.png and
      LABEL_OUT_diff.png
The root mean square error is printed in the standard output.


===============================================================================
====== List of files ==========================================================



The following files implement six image smoothing methods:
    lib_gaussconv.h
    lib_gaussconv.c

The following files contains input/output routines:
    lib_io_gaussconv.h
    lib_io_gaussconv.c

The following files produce binaries for each image
    gaussconv_dct.c
    gaussconv_dft.c
    gaussconv_lindeberg.c
    gaussconv_sampled_kernel.c

The folloling files produce the binaries to examine the consistency
to semigroup property for each of the five methods. 
    demo_gaussconv_dct.c
    demo_gaussconv_dft.c
    demo_gaussconv_lindeberg.c
They are used in the companion demo.

===============================================================================
== Compilation (Linux) ========================================================

To compile

Type

    make

in the directory where the Makefile is located.


