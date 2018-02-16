/*
The "Garfields" code is a Gaussian random field generator, which can also
generate divergenceless and helical fields.
Copyright (C) 2007  Francisco S. Kitaura

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Francisco Kitaura can be reached at kitaura@mpa-garching.mpg.de
*/

/*
F.S. Kitaura would like to thank Torsten Ensslin and Jens Jasche for useful
contributions for the development of this code.
*/
#include <fftw3.h>

#ifndef fftw_real
#define fftw_real double
#endif
#define re(c) ((c)[0])
#define im(c) ((c)[1])

#define fwd true
#define fourier_sp true
#define inv false
#define real_sp false

void FFT ( int dim, int *sz, int factor, bool direction, fftw_complex *in, fftw_complex *out );

void FFTR2CC ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out );

void FFTC2RC ( int dim, int *sz, int factor, fftw_complex *in, fftw_real *out );

void FFTR2C ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out );

void FFTC2R ( int dim, int *sz, int factor, fftw_complex *in, fftw_real *out );

void CONVPREPR2C ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out );

void CONVPREPC2C ( int dim, int *sz, int factor, fftw_complex *in, fftw_complex *out );

void convolution ( int dim, int *sz, int factor, bool space_ina, bool space_inb, bool space_out, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );
