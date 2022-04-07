/*
 * $Header: /Users/garciamillan/Documents/uni/ImperialCL/DNAfolding/HiC_2kb/contactdomains/results22_replicates/gplist-master/DNA_jets_replicates.c$
 * COMPILE: cc -O3 -Wall -o Jets DNA_panoramic.c gplist_util2.c -lm -lgsl -lgslcblas
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

int read_depth(struct data_struct *dt);
int read_hic(const char *filename, struct gplist_header_struct **data, int chr);
int read_expected(const char *basenameEXP, double **data, int chr);
int protractor_grid(double rad, const double resolution, int n, struct point_protractor **point);
int panoramic_curve(struct data_struct *dt, int k, struct point_protractor *point, int num_points, int n, const double resolution);
int panoramic_condition(struct condition_struct *condition, double tolerance_whitelines, int n);

double angle(double x, double y){
    double theta;
    if(y>0) {
     if(x>0) theta=atan(y/x);
     else if(x==0) theta=M_PI_2;
     else theta=M_PI+atan(y/x);}
    else theta= M_PI+atan(y/x);
    return theta;
}

int sector(double x, double y, int n){
    if ((x==0) && (y==0)) return -1;
    double theta;
    theta=angle(x,y);
    if (floor((1.25-theta/M_PI))==1) return n-1;
    return floor(n*(1.25-theta/M_PI)); //this way s increases as x increases
}

int main(int argc, char *argv[]){
    setlinebuf(stdout);
    double rad;
    int n,num_jet_candidates;
    char buffer[8192];
    double x0, x1;
    char chrID[100];
    int k, chr_num, chr_num_read, chr_num_prev;
    int i;
    FILE *fg=stdin;
    int ch;
    int replicate;
    int num_replicates=(int)(sizeof(data)/sizeof(*data));
    int num_conditions=(int)(sizeof(condition)/sizeof(*condition));
    struct point_protractor *point=NULL;

    rad=RADIUS_PROTRACTOR;
    n=NUM_SECTORS;

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
    
    while ((ch = getopt(argc, argv, "r:n:h")) != -1)
            switch (ch) {
                case 'r':
                    if ((rad=atof(optarg))<0) {
                        fprintf(stderr, "radius r must be non-negative.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'n':
                    if ((n=atoi(optarg))<=0) {
                        fprintf(stderr, "number of sectors  must be positive.\n");
                        exit(EXIT_FAILURE);
                    }
                    break;
                default:
                    fprintf(stderr, "Flag %c unknown.\n", ch);
                case 'h':
                    fprintf(stderr, "Usage: ./dna -f filename_compartments -p pattern -M -m min:max -b nbins -r resolution file_contact1 file_contact2 file_contact3 ...\n");
                    fprintf(stderr, "-f filename_compartments specified input file (with compartments). By default the program reads from stdin.\n"
                            "-r radius of the circle sector (default 5e5bp)\n"
                            "-R resolution of the lattice (default 5e3bp)\n"
                            "-n number of sectors we consider in half a circle (default 9)\n"
                            "-c chromosome number (by default, n=1).\n"
                            " ./dna -c 2 -f CD69negDPCTCFKOR1R2_eigenvector_KR100kb.bed -r 500000 -R 5000 -n 9\n");
                    return(0);
                    break;
            }
    
    {//two reads to identify jet candidates from the .bed file
    if ((fg=fopen(bed_file, "rt"))==NULL) {
        fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", bed_file, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    printf("# Info: Reading compartments from bedfile %s\n",bed_file);
        
        chr_num_read=chr_num_prev=0;
        num_jet_candidates=0;
        
    //the first read is to find how many chromosomes and jet candidates there are
        int maxIndexOfChromosomes=0;
    printf("# Info: First read.\n");
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
        if(buffer[0]=='c'){
            if (sscanf(buffer,"chr%s", chrID)) {
                if(chrID[0]>='0' && chrID[0]<='9') {
                    chr_num_read=atoi(chrID);
                    maxIndexOfChromosomes=MAX(maxIndexOfChromosomes, chr_num_read);
                    num_jet_candidates++;
                }
            }
        }
    }

        printf("# Info: Total number of jet candidates %d\n"
               "# Info: Total number of chromosomes %d\n",num_jet_candidates,maxIndexOfChromosomes);
        if(num_jet_candidates>=MAX_INDEX_OF_POSITIONS){
            printf("# Warning: Curtailing num_jet_candidates = %i to %i\n", num_jet_candidates, MAX_INDEX_OF_POSITIONS);
            num_jet_candidates=MAX_INDEX_OF_POSITIONS;
        }
            
        for(replicate=0;replicate<num_replicates;replicate++) {
            data[replicate].num_positions=num_jet_candidates;
            MALLOC(data[replicate].position, data[replicate].num_positions);
        }
    rewind(fg);

    k=0;

    printf("# Info: Second read.\n");
    replicate=0;
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL && k<data[0].num_positions) {
        if(buffer[0]=='c'){
            if (sscanf(buffer,"chr%s %lf %lf",chrID,&x0,&x1)) {
                if(chrID[0]>='0' && chrID[0]<='9'){
                    chr_num_read=atoi(chrID);
                
                    for(replicate=0;replicate<num_replicates;replicate++) {
                        data[replicate].position[k].chromosome=chr_num_read;
                        data[replicate].position[k].start=x0;
                        data[replicate].position[k].end=x1;
                        sprintf(data[replicate].position[k].label, "Pos%d", k+1);
                    }
                    k++;
                }
            }
        }
    }
    fclose(fg);

    printf("# Info: Finished reading jet candidates.\n");

    }

    {//load bamfiles (to find white lines)
        printf("# Info: Reading BAM_depth files to find white lines.\n");
        for(replicate=0;replicate<num_replicates;replicate++) read_depth(&(data[replicate]));
        printf("# Info: Finished reading BAM_depth files.\n");
    }

    printf("# Info: Generating grid of points in semi-circle (protractor).\n");
    int num_points = protractor_grid(rad,resolution,n,&point);
    if(num_points<n) {
        fprintf(stderr,"# Error: The semi-circle grid has not enough points:  %d < %d .\n",num_points,n);
        exit(EXIT_FAILURE);
    }
    printf("# Info: The semi-circle grid has %d points.\n",num_points);
    printf("# Info: The semi-circle has radius %g .\n",rad);
    

    double tolerance_whitelines=TOLWHITELINES*num_points/((double)n);

    {
    printf("# Info: Generating panoramic curves.\n");
    chr_num_prev=0;
    int replicate_prev=-1;
    for(replicate=0;replicate<num_replicates;replicate++) {
        printf("# Info: Working on replicate %s\n", data[replicate].ident);
        for(k=0; k<data[replicate].num_positions; k++){
            printf("# Info: Working on candidate %s, in chromosome %d, replicate %s\n",
                   data[replicate].position[k].label,data[replicate].position[k].chromosome,data[replicate].ident);
            
            chr_num=data[replicate].position[k].chromosome;
            if(chr_num!=chr_num_prev || replicate_prev!=replicate) {
                gplist_free(data[replicate].hic);
                free(data[replicate].expected);
                
                //load HiC data with KR normalisation (observed)
                read_hic(data[replicate].basenameHIC, &(data[replicate].hic), chr_num);
            
                //load HiC expected values
                read_expected(data[replicate].basenameEXP, &(data[replicate].expected), chr_num);
                chr_num_prev=chr_num;
                replicate_prev=replicate;
            }

            panoramic_curve(&(data[replicate]),k,point,num_points,n,resolution);
        }
    }
    }
    
    int j;
    
    {//group replicates by condition, calculate mean and std panoramics
        for (i=0;i<num_conditions;i++) panoramic_condition(&(condition[i]),tolerance_whitelines,n);
    }
    
    double label_angle[n];
    for(i=0;i<n;i++) label_angle[i]=90.*(-1.+(1.+2.*(double)i)/(double)n);

    
        printf("# Info: Printing output.\n");
    
    {//print panoramic curves
        printf("PANORAMIC #label  chromosome  start  end  index  angle");
        for (i=0;i<num_conditions;i++) printf("  %s_mean  %s_std",condition[i].ident, condition[i].ident);
        printf("\n");

        for(k=0; k<data[0].num_positions; k++) for(j=0;j<n;j++) {
            printf("PANORAMIC %s  %d  %.16G  %.16G  %d  %.16G",
                   data[0].position[k].label, data[0].position[k].chromosome, data[0].position[k].start, data[0].position[k].end,  j, label_angle[j]);
            for (i=0;i<num_conditions;i++) printf("  %.16G  %.16G",condition[i].position[k].panoramic_mean[j], condition[i].position[k].panoramic_std[j]);
            printf("\n");
        }
    }
    
    {//print GNUPLOT commands
    char filePath[PATH_MAX],fileName[PATH_MAX];
    if (readlink("/proc/self/fd/1", filePath, sizeof(filePath)-1)!=-1) strncpy(fileName, filePath + 73, PATH_MAX - 73);
                
        printf("#GNUPLOT00 FILE = '%s'\n"
            "#GNUPLOT00 set xrange [-90:90]\n"
            "#GNUPLOT00 a=-200\n"
            "#GNUPLOT00 b=300\n"
            "#GNUPLOT00 set yrange [a:b]\n"
            "#GNUPLOT00 set arrow from 45.,a to 45.,b nohead dt '.'\n"
            "#GNUPLOT00 set arrow from -45.,a to -45.,b nohead dt '.'\n"
            "#GNUPLOT00 set xlabel 'Angle (degrees)'\n"
            "#GNUPLOT00 set ylabel 'Strength'\n"
            "#GNUPLOT00 set key out horizontal bottom\n"
            "#GNUPLOT00 set zeroaxis\n"
            "#GNUPLOT00\n"
               "#GNUPLOT00\n", fileName);
        
        for(k=0; k<data[0].num_positions; k++) {
            printf("#GNUPLOT01_%d set term 'svg'\n"
            "#GNUPLOT01_%d set out '%s_gnu.svg'\n"
            "#GNUPLOT02_%d set term post eps color solid\n"
            "#GNUPLOT02_%d set out '%s_gnu.eps'\n",
                   k+1,k+1,data[0].position[k].label,k+1,k+1,data[0].position[k].label);
            printf("#GNUPLOT%d plot '<grep %s[^0-9] '.FILE.'' u 7:8:9 w e pt 4 lc 'black' t \"Wild-type\","
                   " '<grep %s[^0-9] '.FILE.'' u 7:10:11 w e pt 6 lc 'blue' t \"CTCF ko\","
                   " '<grep %s[^0-9] '.FILE.'' u 7:12:13 w e pt 8 lc 'red' t \"RAD21 CTCF ko\"\n",
                   k,data[0].position[k].label,data[0].position[k].label,data[0].position[k].label);
        }

    }
    
    
printf("# Info: c'est fini.\n");

return 0;
}



int read_depth(struct data_struct *dt)
{
    FILE *fg;
    char *p;
    char buffer[8192];
    int chr;
    int locus,depth,index;
    char chrID[MAX_CHROMOSOME_NAME];

    /* Determine number of chromosomes and each chromosome length */
    
    for (chr=0; chr<=MAX_INDEX_OF_CHROMOSOMES; chr++) {
        dt->BAM[chr].maxIndex=-1;
        dt->BAM[chr].chromosomeName[0]=(char)0;
    }
    if ((fg=fopen(dt->filenameBAM, "rt"))==NULL) {
        fprintf(stderr, "Can't open file \"%s\" (%i::%s)\n", dt->filenameBAM, errno, strerror(errno));
        exit(EXIT_FAILURE);
    } else printf("#Info: reading file \"%s\"\n", dt->filenameBAM);

    dt->maxIndexOfChromosomes=-1;
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                chr=0;
                if (sscanf(buffer,"chr%s %d %d",chrID,&locus,&depth)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr=strtol(chrID, &p, 10);
                    else {chr=-1; continue;}
                    /* Ignoring X,Y chromosomes. */
                    //else if(chrID[0]=='X' || chrID[0]=='x') chr =LAST_CHR-1;
                    //else if(chrID[0]=='Y' || chrID[0]=='y') chr =LAST_CHR;
                    //strings like chr17_random will be ignored
                    if (*p!=0)
                        if (!isspace(*p)) {chr=-1; continue;}
                    
                    index=(int)(locus/resolution-CHROMOSOME_OFFSET);
                    /* At this stage, we keep track only of max of chr and max of index. */
                    dt->maxIndexOfChromosomes=MAX(dt->maxIndexOfChromosomes, chr);
                    if ((chr>=0) && (chr<=MAX_INDEX_OF_CHROMOSOMES)) {
                        if (dt->BAM[chr].maxIndex==-1) strcpy(dt->BAM[chr].chromosomeName, chrID);
                        dt->BAM[chr].maxIndex=MAX(dt->BAM[chr].maxIndex,index);
                    }
                }
            }
        }
    
    if (dt->maxIndexOfChromosomes>=MAX_INDEX_OF_CHROMOSOMES) {
        printf("# Warning: Curtailing dt->maxIndexOfChromosomes=%i to %i\n", dt->maxIndexOfChromosomes, MAX_INDEX_OF_CHROMOSOMES);
        dt->maxIndexOfChromosomes=MAX_INDEX_OF_CHROMOSOMES;
    }
    
    /* Now we know the maximum index of the chromosomes and the maximum indeces of those ... */
    for (chr=0; chr<=dt->maxIndexOfChromosomes; chr++) {
        printf("# Info: Chromosome %i is called [%s] has maximum index %i\n", chr, dt->BAM[chr].chromosomeName, dt->BAM[chr].maxIndex);
        if (dt->BAM[chr].maxIndex>0) {
            MALLOC(dt->BAM[chr].depth, dt->BAM[chr].maxIndex+1);
            for(index=0;index<=dt->BAM[chr].maxIndex;index++) dt->BAM[chr].depth[index]=0.;
        } else {
            dt->BAM[chr].depth=NULL;
            strcat(dt->BAM[chr].chromosomeName, "--EMPTY--");
        }
    }
    
    rewind(fg);
    
#warning "Boilerplate code."
    
    while (fgets(buffer, sizeof(buffer)-1, fg)!=NULL) {
            if(buffer[0]=='c'){
                chr=0;
                if (sscanf(buffer,"chr%s %d %d",chrID,&locus,&depth)) {
                    if(chrID[0]>='0' && chrID[0]<='9') chr=strtol(chrID, &p, 10);
                    else {chr=-1; continue;}
                    /* Ignoring X,Y chromosomes. */
                    //else if(chrID[0]=='X' || chrID[0]=='x') chr =LAST_CHR-1;
                    //else if(chrID[0]=='Y' || chrID[0]=='y') chr =LAST_CHR;
                    //strings like chr17_random will be ignored
                    if (*p!=0)
                        if (!isspace(*p)) {chr=-1; continue;}
                    
                    index=(int)(locus/resolution-CHROMOSOME_OFFSET);
                    /* At this stage, we keep track only of max of chr and max of index. */
                    if ((chr>=0) && (chr<=dt->maxIndexOfChromosomes)) {
                        if ((index>=0) && (index<=dt->BAM[chr].maxIndex)) dt->BAM[chr].depth[index]+=depth;
                    }
                    
                }
            }
        }
    /* End of copied code */
    fclose(fg);

    return(0);
}


int read_hic(const char *basenameHIC, struct gplist_header_struct **data, int chr){
    if ((*data)!=NULL) {
        printf("# Warning: Pointer to data is not NULL. Does this need to be free(3)d?\n");
    }

    char filename[PATH_MAX];
    sprintf(filename, "%s%d_%d_5kb.bin.gz", basenameHIC,chr,chr);
    
    printf("# Info: Reading HiC from file %s\n", filename);
    if (((*data)=gplist_read(filename))==NULL) {
        fprintf(stderr, "# Error: Failed to read %s. (%i::%s)\n", filename, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    return 0;
}

int read_expected(const char *basenameEXP, double **data, int chr){
    FILE *fexpected;
    char buffer[8192];
    int i;
    double expected;
    char filename[PATH_MAX];

    if ((*data)!=NULL) {
        printf("# Warning: Pointer to data is not NULL. Does this need to be free(3)d?\n");
    }
    
    MALLOC((*data),LENGHT_EXP);
    sprintf(filename, "%s%d_%d_5kb.dat", basenameEXP,chr,chr);

    
    if ((fexpected=fopen(filename, "rt"))==NULL) {
        fprintf(stderr, "Cannot open file \"%s\" (%i::%s)\n", filename, errno, strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    printf("# Info: Reading expected from file %s, chromosome %d.\n", filename, chr);
    i=0;
    while(fgets(buffer, sizeof(buffer)-1, fexpected)!=NULL) if(sscanf(buffer,"%lf",&expected)) {
        (*data)[i++]=expected;
            if ((i>=LENGHT_EXP) || (i<0)) {
                    fprintf(stderr, "i=%i out of bounds [0,%i]\n", i, LENGHT_EXP);
                    exit(EXIT_FAILURE);
            }
        }
    fclose(fexpected);
    
    return 0;
}


int protractor_grid(double rad, const double resolution, int n, struct point_protractor **point){
    int N=floor(rad/resolution);
    int i=0;
    int sect;
    double x,y,distance;
    
    MALLOC(*point,2*N*N);
    
    for(y=-floor(N/sqrt(2.)); y<=N; y++){

        for(x=-floor(sqrt(N*N-y*y)); x<=y; x++){
            if((distance=sqrt(x*x+y*y))>N) continue;
            if(distance<0.5) continue; //only checking for distance == 0
            sect=sector(x,y,n);
            if((sect<0) || (sect>=n)) continue;
            if(i>=2*N*N) { fprintf(stderr, "Error: i=%d exceeds bound %d\n", i,2*N*N); exit(EXIT_FAILURE); }

            (*point)[i].x=resolution*x;
            (*point)[i].y=resolution*y;
            (*point)[i].sector=sect;
            (*point)[i].distance=distance;
            i++;
            
        }
    }
    
    return i;
}

int panoramic_curve(struct data_struct *dt, int k, struct point_protractor *point, int num_points, int n, const double resolution){
    int i,sector_number,chr_num,epos;
    double x_bp,y_bp,x_bpr,y_bpr,count;
    
    if (dt->hic==NULL) { fprintf(stderr, "Error: dt->hic==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (dt->expected==NULL) { fprintf(stderr, "Error: dt->expected==NULL unexpected in line %i of %s\n", __LINE__, __FILE__); exit(EXIT_FAILURE); }
    if (point==NULL) { fprintf(stderr, "Pointer 'point' is NULL\n"); exit(EXIT_FAILURE); }
    if ((k>= dt->num_positions) || (k<0)) { fprintf(stderr, "k=%i out of bounds [0,%i]\n", k, dt->num_positions-1); exit(EXIT_FAILURE); }
    
    MALLOC(dt->position[k].panoramic, n);
    MALLOC(dt->position[k].whitelines, n);
    
    for(sector_number=0;sector_number<n;sector_number++) dt->position[k].panoramic[sector_number]=dt->position[k].whitelines[sector_number]=0.; //initialisation

    for(i=0;i<num_points;i++) {
        x_bp=point[i].x+dt->position[k].start;
        y_bp=point[i].y+dt->position[k].end;
        x_bpr=floor(x_bp/resolution);
        y_bpr=floor(y_bp/resolution);
        epos=(int)floor(fabs((x_bp-y_bp)/resolution));

        sector_number=(int)point[i].sector;
        if ((sector_number>=n) || (sector_number<0)) { fprintf(stderr, "sector_number = %i out of bounds [0,%i]\n", sector_number, n); exit(EXIT_FAILURE); }
        if ((epos<0) || (epos>=LENGHT_EXP)) { fprintf(stderr, "Error: epos=%i out of bounds 0..%i\n", epos, LENGHT_EXP);  exit(EXIT_FAILURE); }
        
        count = gplist_fetch_count(dt->hic,x_bpr,y_bpr) - dt->expected[epos];
        dt->position[k].panoramic[sector_number]+= count;

        chr_num=dt->position[k].chromosome;
        
        if (((int)(x_bp/resolution-CHROMOSOME_OFFSET) > dt->BAM[chr_num].maxIndex) || ((int)(x_bp/resolution-CHROMOSOME_OFFSET)<0)) dt->position[k].whitelines[sector_number]++;
        else if(dt->BAM[chr_num].depth[(int)(x_bp/resolution-CHROMOSOME_OFFSET)]==0) dt->position[k].whitelines[sector_number]++;
        
        if (((int)(y_bp/resolution-CHROMOSOME_OFFSET) > dt->BAM[chr_num].maxIndex) || ((int)(y_bp/resolution-CHROMOSOME_OFFSET)<0)) dt->position[k].whitelines[sector_number]++;
        else if(dt->BAM[chr_num].depth[(int)(y_bp/resolution-CHROMOSOME_OFFSET)]==0) dt->position[k].whitelines[sector_number]++;
        
    }
    
    //normalise by total counts
    for(sector_number=0;sector_number<n;sector_number++) dt->position[k].panoramic[sector_number] *= PREF/dt->total_counts;
            
    return 0;
}


int panoramic_condition(struct condition_struct *condition, double tolerance_whitelines, int n){
    int replicate,k,j;
    double m0;
    
    condition->num_positions=condition->data[0]->num_positions;
    MALLOC(condition->position,condition->num_positions);
    
    for(replicate=0; replicate<condition->num_replicates; replicate++){
        condition->data[replicate]->parent= condition;
        printf("#Info %s is daughter of %s\n",condition->data[replicate]->ident,condition->data[replicate]->parent->ident);
    }
    
    for(k=0; k<condition->num_positions; k++){
        MALLOC(condition->position[k].panoramic_mean,n);
        MALLOC(condition->position[k].panoramic_std,n);

        for(j=0;j<n;j++) {
            condition->position[k].panoramic_mean[j] = condition->position[k].panoramic_std[j] = 0.; //initialise
            
            //calculate mean
            m0=0.;
            for(replicate=0; replicate<condition->num_replicates; replicate++){
                if(condition->data[replicate]->position[k].whitelines[j] < tolerance_whitelines){
                    condition->position[k].panoramic_mean[j]+= condition->data[replicate]->position[k].panoramic[j];
                    m0++;
                    }
            }
            if(m0<1) condition->position[k].panoramic_mean[j]= nan("1");
            else condition->position[k].panoramic_mean[j]*= 1./m0;

            //calculate std
            if(m0>1){
                for(replicate=0;replicate<condition->num_replicates;replicate++){
                    if(condition->data[replicate]->position[k].whitelines[j] < tolerance_whitelines) condition->position[k].panoramic_std[j]+= pow(condition->data[replicate]->position[k].panoramic[j] - condition->position[k].panoramic_mean[j],2);
                }
                condition->position[k].panoramic_std[j] = sqrt(condition->position[k].panoramic_std[j]/(m0*(m0-1.)));
            }
    }   }
    
    return 0;
}
