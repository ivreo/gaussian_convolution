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



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3.h>

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

#ifndef M_PI
    #define M_PI 3.14159265359
#endif


#include "lib_gaussconv.h"



////////  Sampled Gaussian kernel - discrete convolution //////////




/** @brief Build mono-dimensional sampled Gaussian kernel
 *   of size 2*rad+1 and standard deviation sigma
 */
double* malloc_gaussian_kernel(double sigma, int rad){
    double* gker = (double*)malloc((2*rad+1)*sizeof(double));
    double tmp = 0;
    
    assert(sigma>=0);
    
    if(sigma>0){
        for(int i = -rad; i <= rad; i++){
            gker[rad+i] = exp(-0.5*(double)i*(double)i/sigma/sigma);
            tmp += gker[rad+i];
        }
        for(int i = -rad; i <= rad; i++)
            gker[rad+i] /= tmp;
    }else{
       for(int i = -rad; i <= rad; i++){
           gker[rad+i] = 0.0;
       }
       gker[rad] = 1.;
    }
    return gker;
}       




/** @brief Apply a convolution with a separable kernel
 * and signal extension by symmetrization
 *
 * @param in             Input image of w X h samples
 * @param out            Output image (same dimension)
 *
 * @param xker           Kernel applied along direction x
 *                          L = r_xker
 *                          w = 2*r_xker+1
 *
 * @param yker           Kernel applied along direction y
 *                          L = r_yker
 *                          w = 2*r_yker+1
 *
 *  Compute:
 *
 *   out[i,j] = \sum_{-r_xker \leq k \leq +r_xker} 
 *                    xker[k]
 *                    . \sum_{-r_xker \leq l \leq +r_xker}
 *                              yker[l].in[i-k,j-l]
 *
 *  Border condition:  symmetrization at border
 *                     (at x = -1./2 and x = w-1./2)
 *
 */
void convolve(double* in, double* out, int w, int h,
              double* xker, int r_xker,
              double* yker, int r_yker){
    int i,j;
    int i_p,j_p;
    int k;
    double tmp;
    
    double* im_tmp = (double*)malloc(w*h*sizeof(double));
                
    /* convolution along x coordinates */
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            tmp=0;
            for(k=-r_xker;k<=r_xker;k++){
                i_p = (i-k+2*h)%(2*h);
                if(i_p>h-1){i_p=2*h-1-i_p;}
              
                tmp += xker[r_xker+k]*in[i_p*w+j];
            }
            im_tmp[i*w+j]=tmp;
        }
    }
    
    /* convolution along y coordinates */
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            tmp=0;
            for(k=-r_yker;k<=r_yker;k++){
                j_p = (j-k+2*w)%(2*w); 
                if(j_p>w-1){j_p=2*w-1-j_p;}
                tmp += yker[r_yker+k]*im_tmp[i*w+j_p];
            }
            out[i*w+j]=tmp;
        }
    }
    
    free(im_tmp);
}


void gaussconv_sampled_kernel(double* in, double* out, int w, int h, double sigma, int k)
{
    int r_gker = (int)(k*sigma);
    double* gker = malloc_gaussian_kernel(sigma, r_gker);
    convolve(in,out,w,h,gker,r_gker,gker,r_gker);
    free(gker);
}



////////  DFT convolution - Gaussian convolution with the DFT interpolation  //

void direct_dft(double* im, fftw_complex* dft, int w, int h)
{
    fftw_plan r2c;
	r2c = fftw_plan_dft_r2c_2d(h, w, im, dft, FFTW_ESTIMATE);
	fftw_execute(r2c);
    double norm = (double)(w*h);
    for(int i=0; i<(w/2+1)*h; i++){
        dft[i][0] /=norm;
        dft[i][1] /=norm;
    }
    fftw_destroy_plan(r2c);
    fftw_cleanup();
}


/** 
 *  dft : size (w/2+1)*h 
 */
void inverse_dft(fftw_complex* dft, double* im, int w, int h)
{
    fftw_plan c2r;
	c2r = fftw_plan_dft_c2r_2d(h, w, dft, im, FFTW_ESTIMATE);
	fftw_execute(c2r);
    fftw_destroy_plan(c2r);
    fftw_cleanup();
}


void gaussian_on_dft(fftw_complex* dftin, fftw_complex* dftout, int w, int h, double sigma)
{
    double s = 2*sigma*sigma*M_PI*M_PI;
	for(int m = 0; m < h; m++){
        int mp = m;
		if (m > -h/2+h-1){mp = m-h;}
        double dmp = (double)mp/(double)h;
		for(int n = 0; n < w/2+1; n++){
			int np = n;
			if (n > -w/2+w-1){ np = n-w;}
            double dnp = (double)np/(double)w;
			double weight = exp(-s*(dmp*dmp + dnp*dnp));
			dftout[m*(w/2+1)+n][0] = dftin[m*(w/2+1)+n][0] * weight;
			dftout[m*(w/2+1)+n][1] = dftin[m*(w/2+1)+n][1] * weight;
			
		}
	}
}


void gaussconv_dft(double* in, double* out, int w, int h, double sigma)
{
	size_t fft_size = (w/2+1)*h;   /// fftw takes advantages of DFT conjugate symmetrie
	fftw_complex* dftin  = (fftw_complex *) fftw_malloc(fft_size * sizeof(fftw_complex));
	fftw_complex* dftout = (fftw_complex *) fftw_malloc(fft_size * sizeof(fftw_complex));
    
    direct_dft(in, dftin, w, h);
    gaussian_on_dft(dftin, dftout, w, h, sigma);
    inverse_dft(dftout, out, w, h);

    fftw_free(dftin);
    fftw_free(dftout);
}





////////  DCT convolution - Gaussian convolution with the DCT interpolation  //



void direct_dct(double* im, double* dct, int w, int h)
{
    fftw_plan r2r;
    r2r = fftw_plan_r2r_2d(h, w, im, dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(r2r);
    double norm = (double)(2*w*2*h);
    for(int i=0; i<w*h; i++)
        dct[i] /=norm;
    fftw_destroy_plan(r2r);
    fftw_cleanup();
}


void inverse_dct(double* dct, double* im, int w, int h)
{
    fftw_plan r2r;
	r2r = fftw_plan_r2r_2d(h, w, dct, im, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(r2r);
    fftw_destroy_plan(r2r);
    fftw_cleanup();
}


void gaussian_on_dct(double* dctin, double* dctout, int w, int h, double sigma)
{
    double s = sigma*sigma*M_PI*M_PI/2.0;
    for(int m=0; m<h; m++){
        double dm = (double)m/(double)h;
        for(int n = 0; n < w; n++){
            double dn = (double)n/(double)w;
            dctout[m*w+n] = dctin[m*w+n]*exp(-s*(dm*dm+dn*dn));
        }
    }
}


void gaussconv_dct(double* in, double* out, int w, int h, double sigma)
{
    double* dctin = (double*)malloc(w*h*sizeof(double));
    double* dctout = (double*)malloc(w*h*sizeof(double));
    direct_dct(in, dctin, w, h);
    gaussian_on_dct(dctin, dctout, w, h, sigma);
    inverse_dct(dctout, out, w, h);
    free(dctin);
    free(dctout);
}







//////////////// Lindeberg's smoothing method

void compute_laplacian_scheme(double* in, double* out, int w, int h, double l)
{   
    assert((0<=l)&&(l <=1));
    int b = 50; // the width of the border

    // extend the image by symmetrization
    int he = h+2*b;
    int we = w+2*b;
    double* iin = (double*)malloc(he*we*sizeof(double));
	for(int i=0;i<he;i++){
		int im = i-b;
		while (im<0) im += 2*h;       
		while (im>2*h-1) im -= 2*h;   
		if (im>h-1) im = 2*h -1 -im;  
		for(int j=0;j<we;j++){
			int jm = j-b;
			while (jm<0){ jm += 2*w;}     
			while (jm>2*w-1){ jm -= 2*w;}
			if (jm>w-1) jm = 2*w -1 -jm; 
			iin[i*we+j] = in[im*w+jm];
		}
	}
    
    // Laplacian
    for(int i=b; i<h+b; i++){
        for(int j=b; j<w+b; j++){
            // \Delta^{+}
            double lapl_p = iin[(i+1)*we+j] + iin[(i-1)*we+j] \
                          + iin[i*we+j+1] + iin[i*we+j-1] \
                          - 4*iin[i*we+j];
            // \Delta^{\times}
            double lapl_t = 0.5*( iin[(i+1)*we+(j+1)] + iin[(i-1)*we+(j+1)] \
                                + iin[(i+1)*we+(j-1)] + iin[(i-1)*we+(j-1)] \
                                - 4* iin[i*we+j]);
            out[(i-b)*w+(j-b)] = l*lapl_p + (1-l)*lapl_t;
        }
    }

    free(iin);
}


void compute_laplacian_scheme_old(double* in, double* out, int w, int h, double l)
{   

    for(int i=0; i <h*w; i++){
        out[i] = 0.0;
    }



    assert((0<=l)&&(l <=1));
    for(int i=1; i<h-1; i++){
        for(int j=1; j<w-1; j++){
            // \Delta^{+}
            double lapl_p = in[(i+1)*w+j] + in[(i-1)*w+j] \
                          + in[i*w+j+1] + in[i*w+j-1] \
                          - 4*in[i*w+j];
            // \Delta^{\times}
            double lapl_t = 0.5*( in[(i+1)*w+(j+1)] + in[(i-1)*w+(j+1)] \
                                + in[(i+1)*w+(j-1)] + in[(i-1)*w+(j-1)] \
                                - 4* in[i*w+j]);
            out[i*w+j] = l*lapl_p + (1-l)*lapl_t;
        }
    }

}


void gaussconv_lindeberg(double* in, double* out, int w, int h, double sigma, double g)
{
    double* lapl = (double*)malloc(w*h*sizeof(double));

    // Euler iteration number
    int P = 4*(int)( (8*(1-g/2)*sigma*sigma) +1);
    // Euler step size
    double dt = sigma*sigma/2/P;

    // initialization
    for(int i=0; i<w*h; i++){
        out[i] = in[i];
    }

    for(int p =0; p<P; p++){
        compute_laplacian_scheme(out, lapl, w, h, g);

        // Euler iteration
        for(int i=0; i<w*h; i++){
            out[i] = out[i] + dt* lapl[i];
        }
    }
    free(lapl);

}
