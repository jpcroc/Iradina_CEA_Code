/*********************************************************************

    Copyright 2019, Christian Borschel

    This file is part of iradina.

    iradina is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iradina is distributed WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with iradina.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************/


/*********************************************************************/
/* This file contains definitions adapted from Corteo (version
   Corteo20090527).
   The index values defined in this file can be used with the scattering
   matrix indexed with 4 mantissa bits in the reduced energy and the
   reduced impact factor.

   Corteo was written by Francois Schiettekatte.
   Corteo was released under the GNU General Public License as
   published by the Free Software Foundation, version 3.
   
   You may obtain the original corteo source code from:
   http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
*/
/*********************************************************************/

/* matrix index calculation parameters */
#define MINE  0.0000019073486328125f  //  minimum reduced energy = 2^(-19)
#define MAXE  2097152.f // maximum reduced energy = 2^21
#define DIME  640    // matrix dimension: (21-(-19))*2^4 (4 mantissa bits: *2^4)
#define BIASE 1728   // (127-19)*2^4: exponent bias correction while performing division by 2^(-4)
#define SHIFTE 19    // keep 4 of the 23 mantissa bits (23-4=19)

#define MINS  0.00000001490116119384765625f // mimimum reduced impact parameter = 2^(-26)
#define MAXS  64.f   // maximum reduced impact parameter = 2^6
#define DIMS  512    // (6-(-26))*2^4
#define BIASS 1616   // (127-26)*2^4: exponent bias correction while performing division by 2^(-16)
#define SHIFTS 19    // keep 4 of the 23 mantissa bits (23-4=19)

#define MIND  16.f         // minimum energy (eV) for stopping power tables = 2^4
#define MAXD  1073741824.f // maximum energy (eV) for stopping power tables = 2^30 ~ 1 GeV
#define DIMD  416          //  (30-4)*2^4
#define BIASD 2096         // (127+4)*2^4: exponent bias correction while performing division by 2^10
#define SHIFTD 19          // keep 4 of the 23 mantissa bits (23-4=19)

