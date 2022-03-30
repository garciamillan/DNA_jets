/*
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results22_replicates/gplist-master/DNA_jets_replicates.c$
 * COMPILE: cc -O3 -Wall -o BAM pseudoBAMfile.c gplist_util2.c -lm -lgsl -lgslcblas
 * EXECUTE: ./Jets > test.dat
 * EXECUTE: ./Jets -r 500000 -R 5000 -n 17 > test.dat
 * date: Fri 17 Apr 2020 13:19:16 BST
 *
 * To generate the binary files, run
 * ./DNAsliding_window -w 5000:5000 -C 17:17 ../RUN04_contact_domains/CD69negDPCTCFKOR1R2/CD69negDPCTCFKOR1R2_ch_17_17_5kb.dat
 * /bin/mv -f gplist_data.bin CD69negDPCTCFKOR1R2_ch_17_17_5kb_coarsened5kb.bin
 * gzip CD69negDPCTCFKOR1R2_ch_17_17_5kb_coarsened5kb.bin
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include "gplist_util2.h"
#include "data_structure.h"

int chr_length[]={
    0,
    249235000,
    243185000,
    197955000,
    191040000,
    180900000,
    171050000,
    159125000,
    146300000,
    141145000,
    135520000,
    134945000,
    133840000,
    115105000,
    107285000,
    102505000,
    90290000,
    81180000,
    78015000,
    59115000,
    62965000,
    48115000,
    51235000
};

int *hic_reads_per_position;
int read_hic(const char *basenameHIC, struct gplist_header_struct **data, int chr);

int main(int argc, char *argv[]){
    setlinebuf(stdout);
    int chr_num;
    int i;
    int ch;
    
    struct gplist_header_struct *data_hic=NULL;



    printf("# $Header$\n");
    printf("# Command:");
    for (i=0; i<argc; i++) printf(" %s", argv[i]);
    printf("\n");
    { time_t tm;
      tm=time(NULL);
      printf("# Info: %s", ctime(&tm));
    }
    
    { char hostname[128];
        gethostname(hostname, sizeof(hostname)-1);
        hostname[sizeof(hostname)-1]=(char)0;
        printf("# Info: Hostname: %s\n", hostname);
    }
    
    { char cwd[1024];
        cwd[0]=0;
        if(getcwd(cwd, sizeof(cwd)-1)!=NULL){
        cwd[sizeof(cwd)-1]=(char)0;
        printf("# Info: Directory: %s\n", cwd);}
    }
    
    printf("# Info: PID: %i\n", (int)getpid());
    

    char basenameHIC[PATH_MAX];
    
    while ((ch = getopt(argc, argv, "s:")) != -1)
            switch (ch) {
                case 's':
                    strncpy(basenameHIC,optarg,PATH_MAX);
                break;
                default:
                    fprintf(stderr, "Flag %c unknown.\n", ch);
            }
    
    

    int num_chr=(int)(sizeof(chr_length)/sizeof(*chr_length));;
    int x_bpr,y_bpr,max_ind_chr;
    double count;

    printf("# Info: Reading HIC data fom bin files.\n");
    for(chr_num=1;chr_num<num_chr;chr_num++) {
        printf("# Info: Working on chromosome %d\n", chr_num);
        
        gplist_free(data_hic);
        free(hic_reads_per_position);

        //load HiC data
        read_hic(basenameHIC, &data_hic, chr_num);
        
        printf("# Info: length of chromosome %d is %g\n", chr_num,data_hic->chr_w_length);

        
        printf("# Info: Generating pseudo BAM file.\n");

        max_ind_chr=ceil(data_hic->chr_w_length/resolution);
        MALLOC(hic_reads_per_position,max_ind_chr);
        for(i=0;i<max_ind_chr;i++) hic_reads_per_position[i]=0;
        
        for(x_bpr=0;x_bpr<max_ind_chr;x_bpr++){
            for(y_bpr=0;y_bpr<max_ind_chr;y_bpr++){
                count = gplist_fetch_count(data_hic,x_bpr,y_bpr);
                hic_reads_per_position[x_bpr]+=count;
            }
            printf("chr%d  %.16G  %d\n",chr_num,x_bpr*resolution,hic_reads_per_position[x_bpr]);
    }
    }
    
    
    
printf("# Info: c'est fini.\n");

return 0;
}

int read_hic(const char *basenameHIC, struct gplist_header_struct **data, int chr){
    if ((*data)!=NULL) {
        printf("# Warning: Pointer to data is not NULL. Does this need to be free(3)d?\n");
    }

    char filename[PATH_MAX];
    sprintf(filename, "%s%d_%d_5kb.bin.gz", basenameHIC,chr,chr);
    
    printf("# Info: Reading HiC from file %s\n", filename);
    if (((*data)=gplist_read(filename))==NULL) {
        fprintf(stderr, "# Error: Failed to read %s . (%i::%s)\n", filename, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    return 0;
}

