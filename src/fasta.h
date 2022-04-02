//
// Created by parthajit on 29/2/20.
//

#ifndef __fasta_H__
#define __fasta_H__

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

struct squence{
    char* data;
    int size;
};

void sequence_create(struct sequence* seq, const char* dat_file_name){
    FILE* fp = fopen(dat_file_name, "r");
    assert(fp != NULL);
    char tmpseq[20000];
    seq->size = 0;
    char line[100];
    int len;
    //cout<<"One\n";
    while(fgets(line, sizeof(line), fp) != NULL) {
    //cout<<"coming inside="<<line<<"END"<<endl;
        if(line[0] == '>') continue;
        len = sprintf(tmpseq+seq->size,"%s", line);
      //  cout<<"Line="<<line<<", len="<<len<<"."<<endl;
        seq->size += len - 1; // To discard \n
    }
    //cout<<"Two\n";
    /* Allocate the memory of the struture and copy */

    seq->data = (char*) malloc((seq->size+2)* sizeof(char));
    //cout<<"SEq size = "<<seq->size<<endl;
    sprintf(seq->data,"%s",tmpseq);
    seq->data[seq->size] = '\0';
    //printf("The val=%s.\n",seq->data);
    //printf("SEQ size = %d\n", seq->size);
    assert(fp != NULL);
    fclose(fp);
}
void sequence_destroy(struct sequence* seq){
    free(seq->data);
}

#endif //  __fasta_H__
