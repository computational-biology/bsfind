/*
 * =====================================================================================
 *
 *       Filename:  bioargs.h
 *
 *    Description:  Parameters for the program
 *
 *        Version:  1.0
 *        Created:  Friday 31 December 2021 04:34:47  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */
 
#ifndef  __bioargs_H__
#define  __bioargs_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct args{
      struct{
	    char basename[128];
	    char ext[10];
	    char path[512];
	    char type;
	    char full_name[512];
      }file;
      struct{
	    char occu;
      }bio;
      struct{
	    char os[10];
      }sys;
};


void args_init(struct args* args);


void args_process_argv(int argc, char* argv[], struct args* args, int file_index[], int* file_count);

#endif   /* ----- #ifndef __bioargs_H__  ----- */

