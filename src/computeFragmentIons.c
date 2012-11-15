/*
 * computeFragmentIons.c R package protViz
 * parts taken from msmsseqass.c deltatectra
 *
 * Copyright 2007, 2008, 2009, 2010
 * Christian Panse <cp@fgcz.ethz.ch>
 *
 * This file is part of deltatectra.
 *
 * deltatectra is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA *
*/

#include <stdio.h>

/*

 TODO

 o give the array of amino acids as argument

*/

void computeFragmentIons(int *n_, char **seq_, double *pim_, double *b_,
			 double *y_)
{

    int i;
    double b, y;
    double C_term;
    double N_term;
    double Electron;
    double Hydrogen;
    double M[26];
    int letter;

    C_term = 17.002740;
    N_term = 1.007825;
    Electron = 0.000549;
    Hydrogen = 1.007825;
    M[0] = 71.037110;
    M[1] = 114.534940;
    M[2] = 160.030649;
    M[3] = 115.026940;
    M[4] = 129.042590;
    M[5] = 147.068410;
    M[6] = 57.021460;
    M[7] = 137.058910;
    M[8] = 113.084060;
    M[9] = 0.000000;
    M[10] = 128.094960;
    M[11] = 113.084060;
    M[12] = 131.040480;
    M[13] = 114.042930;
    M[14] = 0.000000;
    M[15] = 97.052760;
    M[16] = 128.058580;
    M[17] = 156.101110;
    M[18] = 87.032030;
    M[19] = 101.047680;
    M[20] = 150.953630;
    M[21] = 99.068410;
    M[22] = 186.079310;
    M[23] = 111.000000;
    M[24] = 163.063330;
    M[25] = 128.550590;

    b = N_term - Electron;
    y = *pim_;

    for (i = 0; i < *n_; i++) {
	letter = seq_[0][i];
	if (64 < letter && letter < 92) {
	    b += M[letter - 65];
	    b_[i] = b;

	    y_[*n_ - i - 1] = y;
	    y -= M[letter - 65];
	}
    }
}
