/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Friday 25 March 2022 05:16:29  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bioio.h"
#include "bioio.h"
#include "rnabp.h"
#include "helix.h"


 


void gen_pymol(FILE* fp, struct nucbp* rnabp, char* cif, int from, int size, char* seq)
{
     
      int othindex;
      fprintf(fp, "select %s%d%s, ", cif, rnabp[from].cifid, rnabp[from].chain);
      for(int k=0; k<size; ++k){
	    fprintf(fp, "(resi %d and chain %s)", 
			rnabp[from + k].cifid,
			rnabp[from + k].chain);

	    othindex = rnabp[from + k].oth_base_index[0];
	    
	    fprintf(fp, "(resi %d and chain %s)", 
			rnabp[othindex].cifid,
			rnabp[othindex].chain);
      }

      fprintf(fp, "\n");
      fprintf(fp, "as spheres,  %s%d%s\n", cif, rnabp[from].cifid, rnabp[from].chain);
      fprintf(fp, "show cartoon,  %s%d%s\n", cif, rnabp[from].cifid, rnabp[from].chain);

      fprintf(fp, "select %s%d%st, ", cif, rnabp[from].cifid, rnabp[from].chain);
      for(int k=0; k<size; ++k){
	    if(seq[k] == 'T'){

		  othindex = rnabp[from + k].oth_base_index[1];

		  fprintf(fp, "(resi %d and chain %s)", 
			      rnabp[othindex].cifid,
			      rnabp[othindex].chain);

	    }
      }
      fprintf(fp, "\n");
      fprintf(fp, "as spheres, %s%d%st\n", cif, rnabp[from].cifid, rnabp[from].chain);
      fprintf(fp, "show cartoon, %s%d%st\n", cif, rnabp[from].cifid, rnabp[from].chain);
      fprintf(fp, "util.cbay %s%d%st\n", cif, rnabp[from].cifid, rnabp[from].chain);

}

int main(int argc, char* argv[])
{
      char base[32];
      char path[512];
      char ext[512];

      char outfile[512];
      char datfile[512];

      char sep[] = "\n \tNCLBW";
      char* token;

      FILE* fp= stdout;
      

      for(int i=1; i<argc; ++i){
	    fname_split(path, base, ext, argv[i]);
	    if(strcmp(ext, ".dat") != 0){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Invalid extn %s\n", __func__, __FILE__, __LINE__, ext);
		  exit(EXIT_FAILURE);
	    }

	    fname_join(outfile, path, base, ".out");
	    fname_join(datfile, path, base, ".dat");
	    int numres = get_numres(outfile);

	    struct fasta seq;
	    fasta_init(&seq, numres);
	    scanfasta(datfile, &seq);
	    printf("%s\n", seq.data);

	    struct nucbp* rnabp;



	    rnabp = (struct nucbp*) malloc ( numres * sizeof(struct nucbp) );
	    if ( rnabp==NULL ) {
		  fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
		  exit (EXIT_FAILURE);
	    }

	    rnabp_scan_out(rnabp, numres, outfile);
	    
	    struct helix* helix;
	    int hlxcount;
	    helix_init(&helix, numres);
        helix_compute(helix, &hlxcount, rnabp, numres);
        char hlx[512];
        for(int j=0; j<hlxcount; ++j){
            int hsize = helix[j].size;
            int hstart = helix[j].i[0];
           // for(int l=0; l<hsize; ++l){
                strncpy(hlx, seq.data + hstart, hsize);
                hlx[hsize] = '\0';
                if(hsize >=3 ){
			if(strchr(hlx, 'T') != NULL && strchr(hlx,'N') == NULL){
			      fprintf(fp, "DATA  %s %d  %s  %c %d %d %s  ", 
					  base,
					  rnabp[hstart].cifid,
					  rnabp[hstart].chain,
					  rnabp[hstart].ins,
					  hstart+1,
					  hsize,
					  hlx);
			      for(int k=0; k<hsize; ++k){
				    if(hlx[k] == 'T'){
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

			      gen_pymol(fp, rnabp, base, hstart, hsize, hlx);
			}
			
		  }
		 // }
		  }
                
        return 0;

	    int len;
	    int size;
	    token = strtok(seq.data, sep);
	    while( token != NULL ){
		  len = token - seq.data;
		  size = strlen(token);
		  if(size >=3 ){
			if(strchr(token, 'T') != NULL){
			      fprintf(fp, "DATA  %s %d  %s  %c %d %d %s  ", 
					  base,
					  rnabp[len].cifid,
					  rnabp[len].chain,
					  rnabp[len].ins,
					  len+1,
					  size,
					  token);
			      for(int k=0; k<size; ++k){
				    if(token[k] == 'T'){
					  int othloc = rnabp[len+k].oth_base_index[1];
					  fprintf(fp, "    %s:%s-%s",
						      rnabp[len+k].resname,
						      rnabp[othloc].resname,
						      rnabp[len+k].name[1]);
					  if(rnabp[len+k].numbp >2){
						fprintf(fp, "*");
					  }
				    }

			      }
			      fprintf(fp, "\n");

			      gen_pymol(fp, rnabp, base, len, size, token);
			}
			
		  }
		  token = strtok(NULL, sep);
	    }
	    
	    helix_free_all(helix);

	    free ( rnabp );
	    rnabp	= NULL;

	    

	    fasta_free(&seq);
      }

}
