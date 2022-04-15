/*
 * =====================================================================================
 *
 *       Filename:  ligsite.c
 *
 *    Description:  Lingand site binding functions.
 *
 *        Version:  1.0
 *        Created:  Monday 04 April 2022 04:03:20  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include "ligsite.h"

double ligsite_phosp_dist(FILE* fp, struct polymer* poly, struct nucbp* rnabp, char* accn, int hstart, int hsize, char* seq)
{
      fprintf(fp, "\n");
      int res1indx, res2indx;
      double dist2;
      for(int i=0; i<hsize; ++i){
	    res1indx = hstart + i;
	    struct atom* atm1 = residue_get_atom(poly->residues + res1indx, "P");
	    for(int j=0; j<hsize; ++j){
		  res2indx = hstart + j;
		  struct atom* atm2 = residue_get_atom(poly->residues + res2indx, "P");
		  if(atm1 != NULL && atm2 != NULL){

			dist2 = distsqr(atm1->center, atm2->center);
		  }else{
			dist2 = -1.0;
		  }
		  fprintf(fp, "%8.2lf  ", dist2);


	    }
	    fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
      fprintf(fp, "---------------------------------------------------\n");

}
void ligsite_comp(FILE* fp, struct nucbp* rnabp, char* accn, int hstart, int hsize, char* seq)
{
      fprintf(fp, "DATA  %s %d  %s  %c %d %d %s  ",
		  accn,
		  rnabp[hstart].cifid,
		  rnabp[hstart].chain,
		  rnabp[hstart].ins,
		  hstart+1,
		  hsize,
		  seq);
      for(int k=0; k<hsize; ++k){
	    if(seq[k] == 'T'){
		  int othloc = rnabp[hstart+k].oth_base_index[1];
		  fprintf(fp, "    %s:%s-%s",
			      rnabp[hstart+k].resname,
			      rnabp[othloc].resname,
			      rnabp[hstart+k].name[1]);
		  if(rnabp[hstart+k].numbp >2){
			fprintf(fp, "*");
		  }
	    }

      }
      fprintf(fp, "\n");
      fprintf(fp, "---------------------------------------------------\n");
}



void ligsite_pymol(FILE* fp, struct nucbp* rnabp, char* accn, int hstart, int hsize, char* seq)
{

      int othindex;
      fprintf(fp, "select %s%d%s, ", accn, rnabp[hstart].cifid, rnabp[hstart].chain);
      for(int k=0; k<hsize; ++k){
	    fprintf(fp, "(resi %d and chain %s)",
			rnabp[hstart + k].cifid,
			rnabp[hstart + k].chain);

	    othindex = rnabp[hstart + k].oth_base_index[0];

	    fprintf(fp, "(resi %d and chain %s)",
			rnabp[othindex].cifid,
			rnabp[othindex].chain);
      }

      fprintf(fp, "\n");
      fprintf(fp, "as spheres,  %s%d%s\n", accn, rnabp[hstart].cifid, rnabp[hstart].chain);
      fprintf(fp, "show cartoon,  %s%d%s\n", accn, rnabp[hstart].cifid, rnabp[hstart].chain);

      fprintf(fp, "select %s%d%st, ", accn, rnabp[hstart].cifid, rnabp[hstart].chain);
      for(int k=0; k<hsize; ++k){
	    if(seq[k] == 'T'){

		  othindex = rnabp[hstart + k].oth_base_index[1];

		  fprintf(fp, "(resi %d and chain %s)",
			      rnabp[othindex].cifid,
			      rnabp[othindex].chain);

	    }
      }
      fprintf(fp, "\n");
      fprintf(fp, "as spheres, %s%d%st\n", accn, rnabp[hstart].cifid, rnabp[hstart].chain);
      fprintf(fp, "show cartoon, %s%d%st\n", accn, rnabp[hstart].cifid, rnabp[hstart].chain);
      fprintf(fp, "util.cbay %s%d%st\n", accn, rnabp[hstart].cifid, rnabp[hstart].chain);

}
