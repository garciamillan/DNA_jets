/* This code is based on 
 * $ Header: /home/ma/p/pruess/.cvsroot/gplist/gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $
 *
 * It's a stripped-down version of gplist_util, but with the major advantage that no DNA related global vars are
 * used, so that all DNA related data are held in the header. This way, we can handle multiple files.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/errno.h>

#include "gplist_util2.h"

#include "gplist_git_stamps.h"

#define GPLIST_END_OF_LIST (-1)

long long int total_malloced=0LL;
#define MALLOC(a,n) {if ((a=malloc(((long long int)(n))*((long long int)sizeof(*(a)))))==NULL) \
   {fprintf(stderr, "Fatal error (%s::%i) malloc(3)ing %lli bytes for variable %s: %i (%s).\n", \
         __FILE__, __LINE__, (long long int)(((long long int)n)*((long long int)sizeof(*(a)))), #a, errno, strerror(errno)); exit(EXIT_FAILURE);} else { \
    total_malloced+=((long long int)( ((long long int)(n))*((long long int)sizeof(*(a))))); printf("#Info: malloc(3)ed %lli bytes for %s, total: %lli (%g GiB)\n", (long long int)(((long long int)(n))*((long long int)sizeof(*(a)))), #a, total_malloced, ((double)total_malloced)/((double)(1024LL*1024LL*1024LL))); }}


/* gplist stuff from here */


static FILE *general_open(char *filename, int *is_popen);


#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)

double gplist_fetch_count(struct gplist_header_struct *hdr, int source, int target)
{
int k;

if ((k=gplist_fetch_edge(hdr, source, target))!=GPLIST_END_OF_LIST) 
  return(hdr->edge_pool[k].count);
else
  return(0.);
}


/* source and target in units of resolution. */
int gplist_fetch_edge(struct gplist_header_struct *hdr, int source, int target)
{
int first, last, try;

#define CHECK_BOUNDS(a) if (!(GPLIST_WITHIN_BOUNDS_LOCATION(hdr,a))) {return(GPLIST_END_OF_LIST);}

/* I think it's worth the effort to check things... */
CHECK_BOUNDS(source);
CHECK_BOUNDS(target);

#define RETURN(a) {return(a);}

/* Look for the target in the range 
 * gplist_location[source].link_sorted to gplist_location[source].link_sorted+gplist_location[source].num_neighbours-1;
 */
first=hdr->location[source].link;
last=hdr->location[source].link+hdr->location[source].num_neighbours-1;

while (last>=first) {
  try=(first+last)/2;
  //printf(" Trying %i which is %i for target %i\n", try, hdr->edge_pool[try].target, target);
  if (hdr->edge_pool[try].target==target) RETURN(try)
  else if (hdr->edge_pool[try].target>target) last=try-1;
  else first=try+1;
}
RETURN(GPLIST_END_OF_LIST);
#undef RETURN
}


struct gplist_header_struct *gplist_read(char *filename)
{
FILE *file;
struct gplist_header_struct *hdr;
int objects_read;
int is_popen=0;

printf("# Info: Version of git_version_string to follow.\n");
printf("%s", git_version_string);
printf("# Info: CVS Header of this file: $Header: /home/ma/p/pruess/.cvsroot/gplist/gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $\n");
printf("# Info: CVS/git ID of this file: $Id: gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $\n");

if (filename[0]==(char)0) {
  file=stdin;
} else {
  /* if ((file=fopen(filename, "r"))==NULL) */
  if ((file=general_open(filename, &is_popen))==NULL) {
    fprintf(stderr, "# Error: Failed to open %s [%s] for reading, error %i::%s\n", (is_popen) ? "stream" : "file", filename, errno, strerror(errno));
    return(NULL);
  }
}
MALLOC(hdr, 1);
hdr->location=NULL;
hdr->edge_pool=NULL;

/* Could call gplist_free(hdr) to free from here. */

#define FREAD(a,b,c,d,e) if ( (objects_read=fread(a,b,c,d))!=(c) ) { fprintf(stderr, "Expected to read %i objects, read %i objects in line %i of %s. (%i::%s)\n", (c), objects_read, __LINE__, __FILE__, errno, strerror(errno)); {e}; return(NULL); }

FREAD(hdr, sizeof(*hdr), 1, file, free(hdr););

if (strcmp(hdr->magic, GPLIST_MAGIC)) {
  fprintf(stderr, "GPLIST_MAGIC in header not found.\n");
  free(hdr);
  return(NULL);
}

printf("# Info: Header of file [%s]\n", filename);
printf("# Info: Version: %s\n", hdr->version);
printf("# Info: link_num=%i\n", hdr->link_num);
printf("# Info: max_num_neighbours=%i\n", hdr->max_num_neighbours);
printf("# Info: edge_pool_count=%i\n", hdr->edge_pool_count);
printf("# Info: chr_w_length=%g\n", hdr->chr_w_length);
printf("# Info: Delta_w=%g\n", hdr->Delta_w);
printf("# Info: generating_command: %s\n", hdr->generating_command);
printf("# Info: header_size=%i versus %i\n", hdr->header_size, (int)sizeof(struct gplist_header_struct));
printf("# Info: edge_size=%i versus %i\n", hdr->edge_size, (int)sizeof(gplist_edge_strct));
printf("# Info: location_size=%i versus %i\n", hdr->location_size, (int)sizeof(gplist_location_strct));
printf("# Info: spares: %i in gplist_edge_strct and %i in gplist_location_strct\n", (int)sizeof(hdr->location->spare), (int)sizeof(hdr->edge_pool->spare));


if (hdr->link_num!= (((int)(hdr->chr_w_length/hdr->Delta_w))+1) ) {
  fprintf(stderr, "hdr->link_num=%i != ((int)hdr->chr_w_length/hdr->Delta_w)+1)=%i\n", hdr->link_num, ((int)(hdr->chr_w_length/hdr->Delta_w))+1);
  free(hdr);
  return(NULL);
}

MALLOC(hdr->location, hdr->link_num);
MALLOC(hdr->edge_pool, hdr->edge_pool_count);

FREAD(hdr->location, sizeof(*(hdr->location)), hdr->link_num, file, free(hdr->location); free(hdr->edge_pool); free(hdr); );
FREAD(hdr->edge_pool, sizeof(*(hdr->edge_pool)), hdr->edge_pool_count, file, {free(hdr->location); free(hdr->edge_pool); free(hdr);} );

printf("# Info: Read %i snippets with %i edges.\n", hdr->link_num, hdr->edge_pool_count);
/* 
 * Check that you are at the end of the file, unless it's stdin (in which case there may be more data [including headers] or interactive input)
 */
if (file!=stdin) {
  int ch;
  int bytes=0;

/* feof is not triggered until I tried to read data. */
  while ((ch=fgetc(file))!=EOF) bytes++;
  printf("# Info: Testing for trailing data: %i bytes (ch=%i)\n", bytes, ch);
  /* if (!feof(file)) */
  if (bytes!=0) {
    fprintf(stderr, "# Info: Trailing data (%i bytes) unexpected.\n", bytes);
  }
  fclose(file);
}
return(hdr);
}


int gplist_free(struct gplist_header_struct *hdr)
{
if (hdr==0) return(-1);
if (hdr->location) free(hdr->location);
if (hdr->edge_pool) free(hdr->edge_pool);
free(hdr);
return(0);
}



/* As of GP's misc/MannaAdjMoment.c, Header: /home/ma/p/pruess/.cvsroot/misc/MannaAdjMoment.c,v 1.1 2014/06/20 12:29:53 pruess Exp */
#define GZIP_EXTENSION ".gz"
#define BZIP2_EXTENSION ".bz2"
#define BUFFER_SIZE (4096)


#define POPEN(a,b,c) {char command[BUFFER_SIZE]; \
  sprintf(command, "%s %s", b, c); \
  if ((a=popen(command, "r"))==NULL) fprintf(stderr, "# Cannot popen %s for reading. %i::%s\n", command, errno, strerror(errno));}


FILE *general_open(char *filename, int *is_popen)
{
FILE  *file=NULL;

if (strlen(filename)>=strlen(GZIP_EXTENSION))
   if (strcmp(filename+strlen(filename)-strlen(GZIP_EXTENSION), GZIP_EXTENSION)==0) {
     POPEN(file, "gunzip -c", filename);
     if (file!=NULL) {
       *is_popen=1;
       return(file);
     }
   }
if (strlen(filename)>=strlen(BZIP2_EXTENSION))
   if (strcmp(filename+strlen(filename)-strlen(BZIP2_EXTENSION), BZIP2_EXTENSION)==0) {
     POPEN(file, "bunzip2 -c", filename);
     if (file!=NULL) {
       *is_popen=1;
       return(file);
     }
   }
if ((file=fopen(filename, "r"))==NULL) {
  fprintf(stderr, "# Cannot fopen %s for reading. %i::%s\n", filename, errno, strerror(errno));
  }
*is_popen=0;
return(file);
}




