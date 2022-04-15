/*
 * =====================================================================================
 *
 *       Filename:  site.h
 *
 *    Description:  This header is for ligand binding site detection.
 *
 *        Version:  1.0
 *        Created:  Monday 04 April 2022 03:04:00  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#ifndef  __ligsite_H__
#define  __ligsite_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "rnabp.h"
#include "polymer.h"


double ligsite_phosp_dist(FILE* fp, struct polymer* poly, struct nucbp* rnabp, char* accn, int hstart, int hsize, char* seq);
void ligsite_comp(FILE* fp, struct nucbp* rnabp, char* accn, int hstart, int hsize, char* seq);

void ligsite_pymol(FILE* fp, struct nucbp* rnabp, char* accn, int start, int hsize, char* seq);

#endif   /* ----- #ifndef __ligsite_H__  ----- */
