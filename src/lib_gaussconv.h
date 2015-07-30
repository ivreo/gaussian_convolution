#ifndef _LIB_HOWTO_GAUSS_H_
#define _LIB_HOWTO_GAUSS_H_

void gaussconv_sampled_kernel(double* in, double* out, int w, int h, double sigma, int k);

void gaussconv_dft(double* in, double* out, int w, int h, double sigma);

void gaussconv_dct(double* in, double* out, int w, int h, double sigma);

void gaussconv_lindeberg(double* in, double* out, int w, int h, double sigma, double g);

#endif
