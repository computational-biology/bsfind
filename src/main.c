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


#include "biodefs.h"
#include "bioargs.h"
#include "bioio.h"
#include "polymer.h"
#include "rnabp.h"
#include "helix.h"
#include "ligsite.h"






int main(int argc, char* argv[])
{


      char accn[64];

      char outfile[512];
      char datfile[512];
      char sitefile[512];

      char sep[] = "\n \tNCLBW";
      char* token;

      FILE* fp;


      struct args args;

      args_init(&args);

      int* file_index = (int*) malloc(argc * sizeof(int));

      if ( file_index == NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }

      int file_count = 0;
      args_process_argv(argc, argv, &args, file_index, &file_count);
      struct atom* atoms;
      int numatoms;



      for(int i=0; i<file_count; ++i){
	    strcpy(args.file.full_name, argv[file_index[i]]);


	    fname_split(args.file.path, args.file.basename, args.file.ext, args.file.full_name);
	    //	    fname_split(path, base, ext, argv[i]);
	    if(strcmp(args.file.ext, ".cif") != 0){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Invalid extn %s\n", __func__, __FILE__, __LINE__, args.file.ext);
		  exit(EXIT_FAILURE);
	    }

	    fname_join(outfile, args.file.path, args.file.basename, ".out");
	    fname_join(datfile, args.file.path, args.file.basename, ".dat");
	    fname_join(sitefile, args.file.path, args.file.basename, ".site");

	    struct atom* atoms;
	    int atom_sz;
	    if( strcmp(args.file.ext, ".cif") == 0 ){
		  scancif(args.file.full_name, is_modi_nucleic, NULL, NULL, &atoms, &atom_sz, NUC_TYPE, NULL, 'S');
	    }else{    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... invalid extension of file. (%s)\n", __func__, args.file.ext);
		  exit(EXIT_FAILURE);
	    }


	    struct polymer polymer;
	    polymer_create(&polymer, atoms, atom_sz);

	    




	    fp	= fopen(sitefile, "w" );
	    if ( fp == NULL ) {
		  fprintf ( stderr, "couldn't open file '%s'; %s\n",
			      sitefile, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }
	    fprintf(fp, "===================================================\n");
	    int numres = get_numres(outfile);

	    struct fasta seq;
	    fasta_init(&seq, numres);
	    scanfasta(datfile, &seq);
	    //	    printf("%s\n", seq.data);

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
			if(strchr(hlx, 'T') != NULL && strchr(hlx,'N') == NULL && strchr(hlx, 'W') == NULL){
			      strcpy(accn, args.file.basename);
			      ligsite_comp( fp, rnabp, accn, hstart, hsize, hlx);
			      ligsite_phosp_dist(fp, &polymer, rnabp, accn, hstart, hsize, hlx);
			      ligsite_pymol(fp, rnabp, accn, hstart, hsize, hlx);
			      fprintf(fp, "===================================================\n");
			}

		  }
	    }

	    fasta_free(&seq);

	    helix_free(helix);

	    free ( rnabp );
	    rnabp	= NULL;



	    if( fclose(fp) == EOF ) {			/* close output file   */
		  fprintf ( stderr, "couldn't close file '%s'; %s\n",
			      sitefile, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }
	    free(atoms);
	    atoms = NULL;
	    polymer_free(&polymer);
      }

}
