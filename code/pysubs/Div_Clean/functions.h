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

/* */

int factorize ( int dim, int *sz );

double min ( int factor, double *in );

double max ( int factor, double *in );

fftw_real factorial(int term);

void complex_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void complexfactor_mult ( int factor, fftw_real in_a, fftw_complex *in_b, fftw_complex *out );

void complex_cfactor_mult ( int factor, fftw_complex in_a, fftw_complex *in_b, fftw_complex *out );

void complex_mult1 (  fftw_complex in_a, fftw_complex in_b, fftw_complex out );

void realcomplex_mult ( int factor, fftw_real *in_a, fftw_complex *in_b, fftw_complex *out );

void complex_vec_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex out );

void complexreal_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void real_zero ( int factor, fftw_real *out );

void real_mult ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out );

void real_identity ( int factor, fftw_real *out );

void complex_identity ( int factor, fftw_complex in, fftw_complex *out );

void complex_zeropad ( int *sz, int factor, fftw_complex *in, fftw_complex *out );

void complex_zero ( int factor, fftw_complex *out );

void complex_eq ( int factor, fftw_complex *in, fftw_complex *out );

void complex_eq1 ( fftw_complex in, fftw_complex out );

fftw_real real_average ( int factor, fftw_real *in );

void real_delta ( int factor, fftw_real *in, fftw_real *out );

void realcomplex_eq ( int factor, fftw_real *in, fftw_complex *out );

void realcomplex1_eq ( fftw_real in, fftw_complex out );

void complexreal_eq ( int factor, fftw_complex *in, fftw_real *out  );

void complexreal1_eq ( fftw_complex in, fftw_real out  );

void complex_sum ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void complex_subtract ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void complexreal_sum ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void real_sum ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out );

void real_subtract ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out );

void real_eq ( int factor, fftw_real *in, fftw_real *out );

void realfactor_mult ( int factor, fftw_real in_a, fftw_real *in_b, fftw_real *out );

fftw_real  action (int factor, fftw_real *in_a,  fftw_real *in_b);

void power_spectrum ( int factor, fftw_complex *in, fftw_real *out );

void power_fourier ( int dim, int *sz, int factor, fftw_real *in_s, fftw_real *out );

void inverse ( int factor, fftw_real *in, fftw_real *out );

void rootdiag ( int factor,  fftw_real *in, fftw_real *out );

void complex_inverse (  int factor, fftw_complex *in, fftw_complex *out );

void complex_inverse1 ( fftw_complex in, fftw_complex out );

void conjugate ( int factor, fftw_complex *in, fftw_complex *out );

void conjugate1 ( fftw_complex in, fftw_complex out );

void FFT ( int dim, int *sz, int factor, bool direction, fftw_complex *in, fftw_complex *out );

void FFTR2C ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out );

void convolution ( int dim, int *sz, int factor, bool space_ina, bool space_inb, bool space_out, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out );

void hermRedun ( int dim, int *sz, int factor, fftw_complex *in, fftw_complex *out);

fftw_real  variance ( int factor, fftw_real * in, fftw_real mu);

void HERMITICITY(fftw_complex *A_C,int N1, int N2, int N3);

