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

#include <string>
using namespace std;

extern	bool clean_div;
extern	int  Nx, Ny, Nz, seed;
extern	double Lx, Ly, Lz, Vol;
extern	double alpha, PowIndex, kmin, kmax, Bvar;
extern	string  out_file; 