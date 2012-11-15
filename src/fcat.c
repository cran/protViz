#include <stdio.h>

/*
 *  fcat.c
 *
 * Copyright 2006
 * Christian Panse <cp@fgcz.ethz.ch>
 *
 * This file is part of the R-package protViz.
 *
 * protViz is free software; you can redistribute it and/or modify
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




/*
 * Christian Panse <cp@fgcz.ethz.ch> 
 * Fri Sep 29 14:35:42 CEST 2006
 *
 * DESCRIPTION:
 * The program counts the occurrence of characters of the input stream.
 * 
 * COMPILER OPTION:
 * gcc -o fcat fcat.c -O3 -pedantic -Wall
 *
 * EXAMPLE:
 * ./fcat < fgcz_swissprot_20041208.fasta > /tmp/o
 *
 */


int fcat()
{

    int i;

    while ((i = getchar()) != EOF) {
	if (i == 62) {
	    while ((i = getchar()) != '\n');
	    printf("\n");
	}
	else if (i == '\n') {}
	else{
    	   printf("%c", i);
	}
    }
   printf("\n");


    return 0;
}
