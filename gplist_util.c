/*
 * $Header: /home/ma/p/pruess/.cvsroot/gplist/gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/errno.h>

#include "gplist_util.h"

#define GPLIST_END_OF_LIST (-1)

struct gplist_header_struct gplist_header;
double gplist_chr_w_length;
double gplist_Delta_w;

long long int total_malloced=0LL;
#define MALLOC(a,n) {if ((a=malloc(((long long int)(n))*((long long int)sizeof(*(a)))))==NULL) \
   {fprintf(stderr, "Fatal error (%s::%i) malloc(3)ing %lli bytes for variable %s: %i (%s).\n", \
         __FILE__, __LINE__, (long long int)(((long long int)n)*((long long int)sizeof(*(a)))), #a, errno, strerror(errno)); exit(EXIT_FAILURE);} else { \
    total_malloced+=((long long int)( ((long long int)(n))*((long long int)sizeof(*(a))))); printf("#Info: malloc(3)ed %lli bytes for %s, total: %lli (%g GiB)\n", (long long int)(((long long int)(n))*((long long int)sizeof(*(a)))), #a, total_malloced, ((double)total_malloced)/((double)(1024LL*1024LL*1024LL))); }}

#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) ) 
#define MIN(a,b) ( ((a)<(b)) ? (a) : (b) ) 

/* gplist stuff from here */


gplist_location_strct *gplist_location;
int *gplist_link; 
int **gplist_link_last;
gplist_edge_strct *gplist_edge_pool_sequential;
gplist_edge_strct *gplist_edge_pool_sorted;
int *gplist_visited;
int gplist_link_num;
int gplist_edge_pool_sequential_num;
int gplist_edge_pool_sequential_count;
int gplist_edge_pool_sorted_count;
int gplist_max_num_neighbours;

#define GPLIST_CORR_WINDOW_SIZE (20)
double *gplist_window_memory[GPLIST_CORR_WINDOW_SIZE];
int gplist_window_memory_pos;
double *gplist_window_sum;


int gplist_num_newlines(char *filename, long long int *max_bp);


int gplist_init(double chr_w_length, double Delta_w, char *filename)
{
int source;
int max;
long long int num_newlines;
long long int max_bp;

num_newlines=gplist_num_newlines(filename, &max_bp);
if (chr_w_length>0) {
  if (max_bp>chr_w_length) {
    fprintf(stderr, "# Error: largest coordinate %lli found in %s exceed length supplied %g.\n", max_bp, filename, chr_w_length);
    exit(EXIT_FAILURE);
  }
} else chr_w_length=max_bp;

max=chr_w_length/Delta_w;
gplist_link_num=max+1;



gplist_chr_w_length=chr_w_length;
gplist_Delta_w=Delta_w;
/* malloc gplist_edge_pool_sequential, gplist_link and initialise gplist_link */ 
MALLOC(gplist_link, gplist_link_num);
MALLOC(gplist_visited, gplist_link_num);
//gplist_edge_pool_sequential_num=gplist_link_num*GPLIST_AVERAGE_EDGES_PER_NODE;
gplist_edge_pool_sequential_num=2*(num_newlines+gplist_link_num);
MALLOC(gplist_edge_pool_sequential, gplist_edge_pool_sequential_num);

/* Only if there is no coarse graining we don't need to scan through our entire list to see whether a 
 * contact count has to be added or can simply be appended. 
 * Hm. Maybe I always do "no_coarse_graining" and merge identical entries afterwards? */

int no_coarse_graining=1;

if (no_coarse_graining) {
  gplist_no_coarse_graining=1;
  MALLOC(gplist_link_last, gplist_link_num);
  for (source=0; source<gplist_link_num; source++) {
    gplist_link_last[source]=&(gplist_link[source]);
  }
}


/* Edges are initialised with each gplist_link having two edges to its nn. */
for (source=0; source<gplist_link_num; source++) {
  gplist_link[source]=GPLIST_END_OF_LIST;
  gplist_visited[source]=0;
  if (GPLIST_NN_LINK_COUNT>0) {
    gplist_add_edge(source, source-1, GPLIST_NN_LINK_COUNT);
    gplist_add_edge(source, source+1, GPLIST_NN_LINK_COUNT);
  }
}


return(0);
}



#define DD fprintf(stderr, "%s::%i\n", __FILE__, __LINE__)


/* view says 18865774 
 * Code in RUN18 says 18865774 
 * New code says says 18865774 */
int gplist_num_newlines(char *filename, long long int *max_bp) 
{
FILE *in;
int c, newlines=0;
//long long int bp1, bp2;
double bp1, bp2;

printf("# Info: File %s has ... ", filename);
fflush(stdout);

if (max_bp!=NULL) *max_bp=0LL;

if ((in=fopen(filename, "rt"))==NULL) {
  fprintf(stderr, "Cannot open %s for reading.\n", filename);
  exit(EXIT_FAILURE);
}
if (max_bp==NULL) {
while ((c=fgetc(in))!=EOF)
  if (c=='\n') newlines++;
} else {
  while (1) {
    //if (fscanf(in, "%Li %Li", &bp1, &bp2)==2) *max_bp=MAX(MAX(bp1,bp2),*max_bp);
    if (fscanf(in, "%lf %lf", &bp1, &bp2)==2) *max_bp=(long long int)MAX(MAX(bp1,bp2),*max_bp);
    while ((c=fgetc(in))!=EOF) {if (c=='\n') {newlines++; break;}}
    if (c==EOF) break;
  }
}

fclose(in);
printf("%i lines.\n", newlines);
return(newlines);
}



int gplist_add_edge(int source, int target, double count)
{
int *p;


/* It's fine to allow that! */
//if (source==target) return(0);
if (count==0.0) return(0);

if (!GPLIST_WITHIN_BOUNDS_LOCATION(source)) return(-1);
if (!GPLIST_WITHIN_BOUNDS_LOCATION(target)) return(-1);

if (gplist_no_coarse_graining) {
  int *last;

  last=gplist_link_last[source];
  *last=gplist_edge_pool_sequential_count++;

  if (gplist_edge_pool_sequential_count>=gplist_edge_pool_sequential_num) {
    fprintf(stderr, "gplist_edge_pool_sequential about to get exhausted. Bailing out. (%i > = %i)\n", gplist_edge_pool_sequential_count, gplist_edge_pool_sequential_num);
    exit(EXIT_FAILURE);
  }
  /* The address of the last next in the linked list is in gplist_link_last. */
  gplist_edge_pool_sequential[*last].next=GPLIST_END_OF_LIST;
  gplist_link_last[source]=&(gplist_edge_pool_sequential[*last].next);
  gplist_edge_pool_sequential[*last].target=target;
  gplist_edge_pool_sequential[*last].count=count;
  gplist_edge_pool_sequential[*last].visited=0;
  gplist_edge_pool_sequential[*last].flags=0;
  return(0);
}

//printf("%i to %i count %g at gplist_edge_pool_sequential_num=%i of %i\n", source, target, count, gplist_edge_pool_sequential_count, gplist_edge_pool_sequential_num);

if (gplist_edge_pool_sequential_count>=gplist_edge_pool_sequential_num) {
  fprintf(stderr, "gplist_edge_pool_sequential about to get exhausted. Bailing out. (%i > = %i)\n", gplist_edge_pool_sequential_count, gplist_edge_pool_sequential_num);
  exit(EXIT_FAILURE);
}


/* Edges are sorted by decreasing count.  --- No, they are not. I wouldn't know how to maintain that. I am updating them (change their count), plus I do not know what order they are arriving at.
 * *p is a pointer to the current entry. 
 * gplist_link is from the gplist_edge_pool_sequential. */

/* Can I maintain the _same_ edge_pool entry for neighbours of different locations? Na. Each has only one decendant. */

/* First find out whether the target exists already. Need to sort afterwards. */

for (p=&(gplist_link[source]); ((*p)!=GPLIST_END_OF_LIST); p=&(gplist_edge_pool_sequential[*p].next)) {
  if (gplist_edge_pool_sequential[*p].target==target) {
    gplist_edge_pool_sequential[*p].count+=count;
    return(0);
  }
}


/* p should point to me and I should point to what p used to point to -- which is bound to be NULL. */
gplist_edge_pool_sequential[gplist_edge_pool_sequential_count].next=*p;
*p=gplist_edge_pool_sequential_count++;
gplist_edge_pool_sequential[*p].target=target;
gplist_edge_pool_sequential[*p].count=count;
gplist_edge_pool_sequential[*p].visited=0;
gplist_edge_pool_sequential[*p].flags=0;

/*
if (gplist_edge_pool_sequential_count==100000) {
  gplist_print_edges();
  exit(0);
}
*/

return(0);





/* First attempt 
 * for (p=gplist_link[source]; ; p=p->next) {
 *   if (p==NULL) {
 *     / * End of list. --- next line broken already. Need to tell parent! * /
 *     p=gplist_edge_pool+gplist_edge_pool_count++;
 *     p->target=gplist_link[target];
 *     p->count=count;
 *     p->next=NULL;
 *     break;
 *   }
 *   if (p->target==target) {
 *     p->count+=count;
 *     break;
 *   }
 *   if (p->count<count) {
 *   XXX Now I need to know the parent.
 * 
 *   }
 * }
 * 
 *  */


}

int gplist_compare_target(const void *av, const void *bv)
{


return(((gplist_edge_strct *)av)->target - ((gplist_edge_strct *)bv)->target);
/*
gplist_edge_strct *a, *b;

a=(gplist_edge_strct *)av;
b=(gplist_edge_strct *)bv;

return(a->target - b->target);
*/
}


int gplist_reorganise(void)
{
int source;
gplist_edge_strct *sort_pool;
int p, i;
int merges;
int *latest_next;

/* organise the data into gplist_location_strct */
MALLOC(gplist_location, gplist_link_num);


gplist_max_num_neighbours=-1;
for (source=0; source<gplist_link_num; source++) {
  gplist_location[source].link_sequential=gplist_link[source];
  gplist_location[source].visited=gplist_visited[source];
  gplist_location[source].num_neighbours=0;
  gplist_location[source].total_count=0.;
  gplist_location[source].flags=0.;
  gplist_location[source].corr=0.;

  for (gplist_location[source].num_neighbours=0, p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next) {
    gplist_location[source].num_neighbours++;
    gplist_location[source].total_count+=gplist_edge_pool_sequential[p].count;
  }
  if (gplist_max_num_neighbours<gplist_location[source].num_neighbours)
    gplist_max_num_neighbours=gplist_location[source].num_neighbours;
}
printf("# Info: gplist_max_num_neighbours=%i\n", gplist_max_num_neighbours);

/* Sort by target */
MALLOC(sort_pool, gplist_max_num_neighbours);

printf("# Info: Sorting neighbours of %i snippets\n", gplist_link_num);

for (source=0; source<gplist_link_num; source++) {
  if (gplist_location[source].num_neighbours==0) continue;

  merges=0;

  for (i=0, p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next, i++) {
    sort_pool[i]=gplist_edge_pool_sequential[p];
  }
  if (i!=gplist_location[source].num_neighbours) {
    fprintf(stderr, "ERROR: found %i targets, expected %i\n", i, gplist_location[source].num_neighbours);
    exit(EXIT_FAILURE);
  }
  /* I don't know why the code was getting stuck. It felt like it was in an infinite loop. */
  qsort(sort_pool, gplist_location[source].num_neighbours, sizeof(*sort_pool), gplist_compare_target);
  if ((source % 100)==0)
    fprintf(stdout, " # Info: %i of %i sorted.\n", source, gplist_link_num);

  /* Precaution in case of gplist_no_coarse_graining -- check for repeated entries. 
   * Begin at the end. Say the targets are:
   * ABBCCCD
   * You compare the last with the one before last, copy the counts over if necessary.
   * Destroy one of them by adjusting the next pointers.
   * Then re-order.
   *
   * */
  if (gplist_no_coarse_graining) {

    for (i=0; i<gplist_location[source].num_neighbours-1; i++)
      sort_pool[i].next=i+1;
    sort_pool[i].next=GPLIST_END_OF_LIST;

    for (i=gplist_location[source].num_neighbours-1; i>0; i--) {
      if (sort_pool[i].target==sort_pool[i-1].target) {
	/* Merge them */
	sort_pool[i-1].count+=sort_pool[i].count;
	sort_pool[i-1].next=sort_pool[i].next;
	merges++;
      }
    }
    if (merges) 
    {
      gplist_location[source].num_neighbours=0;
      for (p=0; p!=GPLIST_END_OF_LIST; p=sort_pool[p].next) {
	//printf("# Info: %i copied to %i\n", p, gplist_location[source].num_neighbours);
	sort_pool[gplist_location[source].num_neighbours++]=sort_pool[p];
      }
      printf("# Info: source %i had %i merges\n", source, merges);
    }
  }


  //#warning "XXX Still need to implement the merging that has happened above!"
  /* Write the result back into the original locations. They are spread all over the place. I need to overwrite target, count, visited, NOT next. */
  /* OLD 
  for (i=0, p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next, i++) {
    gplist_edge_pool_sequential[p].target=sort_pool[i].target;
    gplist_edge_pool_sequential[p].count=sort_pool[i].count;
    gplist_edge_pool_sequential[p].visited=sort_pool[i].visited;
    gplist_edge_pool_sequential[p].flags=sort_pool[i].flags;
  }
  */
  
  /*
  for (i=0, p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next, i++);
  if ((i!=gplist_location[source].num_neighbours) && (merges==0)) {
    fprintf(stderr, "# Error: num_neighbours=%i i=%i\n", gplist_location[source].num_neighbours, i);
  }
*/


  /* Now I copy the gplist_location[source].num_neighbours entries from sort_pool to gplist_edge_pool_sequential
   * and I terminate gplist_edge_pool_sequential with GPLIST_END_OF_LIST. */
  latest_next=NULL;
  for (i=0, p=gplist_location[source].link_sequential; i<gplist_location[source].num_neighbours; p=gplist_edge_pool_sequential[p].next, i++) {
    gplist_edge_pool_sequential[p].target=sort_pool[i].target;
    gplist_edge_pool_sequential[p].count=sort_pool[i].count;
    gplist_edge_pool_sequential[p].visited=sort_pool[i].visited;
    gplist_edge_pool_sequential[p].flags=sort_pool[i].flags;
    latest_next=&(gplist_edge_pool_sequential[p].next);
  }
  if ((p!=GPLIST_END_OF_LIST) && (merges==0)) {
    fprintf(stderr, "# Error: Last gplist_edge_pool_sequential.next is not GPLIST_END_OF_LIST but there were no merges.\n");
    exit(EXIT_FAILURE);
  }
  /* 19 Nov 2019
   * In response to a SIGSEGV that Rosalba got these days:
   * The if clause below didn't used to be here (it was straight 
   * gplist_edge_pool_sequential[p].next=GPLIST_END_OF_LIST). I was terminating
   * with GPLIST_END_OF_LIST irrespective of whether p was a valid index or
   * not. That makes no sense. Only if p is still valid I should do the termination,
   * which can only ever occur if there had been merging. Otherwise the last p is 
   * GPLIST_END_OF_LIST.
   * I can only hope I did not spoil any of the data
   * before.
   * Ah! I think it was gplist_visited that had been malloc before, so it's probably
   * that one that is normally affected by it. 
   * Hang on. No. The last next needs to be GPLIST_END_OF_LIST. If p is valid, I might 
   * write GPLIST_END_OF_LIST there as well, but that's a bit pointless. I need to
   * terminate the list not there, but at the last one.
   *
   *
   * I should add that I have run the code with and without latest_next (only with the if (p!=GPLIST_END_OF_LIST) to avoid
   * the SIGSEGV tripping) and found 
   * no difference in the output for the offending file:
   *  pruess@macomp001:/home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_20181026/RUN18_jettiness> ../DNAsliding_window -w 5000:5000 -C 15:15 /home/clustor2/ma/g/gunsim/DNAfolding/HICdata_2kb_2181026/RUN08_Ccode/CD69negDPDKOR1R2KR/CD69negDPDKOR1R2KR_ch_15_15_5kb.dat
   *
   *  diff gplist_data.bin gplist_data.bin.old01
   */


  
  if (p!=GPLIST_END_OF_LIST) gplist_edge_pool_sequential[p].next=GPLIST_END_OF_LIST;
  if (latest_next!=NULL) {
    if ((*latest_next!=GPLIST_END_OF_LIST) && (merges==0)) {
      fprintf(stderr, "# Error: Last gplist_edge_pool_sequential.next pointed to by latest_next is not GPLIST_END_OF_LIST but there were no merges.\n");
      exit(EXIT_FAILURE);
    }    
    *latest_next=GPLIST_END_OF_LIST;
  }



}

printf("# Info: Done with sorting.\n");
MALLOC(gplist_edge_pool_sorted, gplist_edge_pool_sequential_count);
gplist_edge_pool_sorted_count=0;
for (source=0; source<gplist_link_num; source++) {
  if ((source % 100)==0)
    fprintf(stdout, "# Info: %i of %i rewritten.\n", source, gplist_link_num);
  if (gplist_location[source].num_neighbours==0) {
    gplist_location[source].link_sorted=GPLIST_END_OF_LIST;
    memset(gplist_location[source].spare, 0, sizeof(gplist_location[source].spare));
    continue;
  }
  gplist_location[source].link_sorted=gplist_edge_pool_sorted_count;
  for (p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next, gplist_edge_pool_sorted_count++) {
    gplist_edge_pool_sorted[gplist_edge_pool_sorted_count]=gplist_edge_pool_sequential[p];
    memset(gplist_edge_pool_sorted[gplist_edge_pool_sorted_count].spare, 0, sizeof(gplist_edge_pool_sorted[gplist_edge_pool_sorted_count].spare));
    if (gplist_edge_pool_sequential[p].next!=GPLIST_END_OF_LIST) gplist_edge_pool_sorted[gplist_edge_pool_sorted_count].next=gplist_edge_pool_sorted_count+1;
    else gplist_edge_pool_sorted[gplist_edge_pool_sorted_count].next=GPLIST_END_OF_LIST;
  }
}


free(sort_pool);
return(0);
}



/* Returns the index within the edge pool of the edge pointing from source to target. 
 * This function relies on the linked list being sorted. */
int gplist_fetch_edge_sequential(int source, int target)
{
int p;


//printf("gplist_fetch_edge_sequential: %i %i", source, target);
//#define RETURN(a) {printf(" %i\n", (a)); return(a);}
#define RETURN(a) {return(a);}

for (p=gplist_location[source].link_sequential; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sequential[p].next) {
  if (gplist_edge_pool_sequential[p].target==target) RETURN(p);
  if (gplist_edge_pool_sequential[p].target>target) RETURN(GPLIST_END_OF_LIST);
}
RETURN(GPLIST_END_OF_LIST);
#undef RETURN
}

/* Ha! This is hilarous. I wrote this function twice... (Second time around, I was wondering whether to put the last return in an else or not...)*/
double XXXgplist_fetch_count_sorted(int source, int target)
{
int edge;

if ((edge=gplist_fetch_edge_sorted(source, target))!=GPLIST_END_OF_LIST) return(gplist_edge_pool_sorted[edge].count);
return(0.);
}

/* XXX Next: Rewrite reorganise to put all edges down consecutively, so that you can use divide and conquer for gplist_fetch_edge. */

int gplist_fetch_edge_sorted(int source, int target)
{
int first, last, try;

#define CHECK_BOUNDS(a) if (!(GPLIST_WITHIN_BOUNDS_LOCATION(a))) {return(GPLIST_END_OF_LIST);}

/* I think it's worth the effort to check things... */
CHECK_BOUNDS(source);
CHECK_BOUNDS(target);

//printf("gplist_fetch_edge_sorted: %i %i has %i neighbours", source, target, gplist_location[source].num_neighbours);
/*
for (p=gplist_location[source].link_sorted; p!=GPLIST_END_OF_LIST; p=gplist_edge_pool_sorted[p].next) {
  printf(" (%i)->%i", p, gplist_edge_pool_sorted[p].target);
}
*/

//#define RETURN(a) {printf(" ... returns %i\n", (a)); return(a);}
#define RETURN(a) {return(a);}

/* Look for the target in the range 
 * gplist_location[source].link_sorted to gplist_location[source].link_sorted+gplist_location[source].num_neighbours-1;
 */
first=gplist_location[source].link_sorted;
last=gplist_location[source].link_sorted+gplist_location[source].num_neighbours-1;

while (last>=first) {
  try=(first+last)/2;
  //printf(" Trying %i which is %i for target %i\n", try, gplist_edge_pool_sorted[try].target, target);
  if (gplist_edge_pool_sorted[try].target==target) RETURN(try)
  else if (gplist_edge_pool_sorted[try].target>target) last=try-1;
  else first=try+1;
}
RETURN(GPLIST_END_OF_LIST);
#undef RETURN
}

double gplist_fetch_count_sorted(int source, int target)
{
int k;

if ((k=gplist_fetch_edge_sorted(source, target))!=GPLIST_END_OF_LIST) 
  return(gplist_edge_pool_sorted[k].count);
else
  return(0.);
}

int gplist_gyration(void) 
{
int i, max_dist, delta;
double m0, m1, m2, count;

for (i=0; i<gplist_link_num; i++) {
  max_dist=MIN(i, gplist_link_num-i-1);
  m0=0.0;
  m1=0.0;
  m2=0.0;
  for (delta=0; delta<=max_dist; delta++) {
    count=gplist_fetch_count_sorted(i-delta,i+delta);
    m0+=count;
    m1+=(count*delta);
    m2+=(count*delta*delta);
  }
#define NORM(a) (((a)!=0.) ? (a) : (-1.) )
  printf("# GYRATION %i %10.20g %10.20g %10.20g\n", i, m0, m1/NORM(m0), m2/NORM(m0));
  for (delta=0; delta<=max_dist; delta++) {
    count=gplist_fetch_count_sorted(i-delta,i+delta);
    printf("# GYRPROFILE %i %i %10.20g %10.20g %10.20g\n", i, delta, count/NORM(m0), count, m0);
  }
}

return(0);
}



int gplist_spears(void) 
{
int i, max_dist, delta;
double m0, m1, m2, count;

#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) ) 
#define MIN(a,b) ( ((a)<(b)) ? (a) : (b) ) 
for (i=0; i<gplist_link_num; i++) {
  max_dist=MIN(i, gplist_link_num-i-1);
  m0=0.0;
  m1=0.0;
  m2=0.0;
  for (delta=0; delta<=max_dist; delta++) {
    count=gplist_fetch_count_sorted(i-delta,i+delta);
    m0+=count;
    m1+=(count*delta);
    m2+=(count*delta*delta);
  }
#define NORM(a) (((a)!=0.) ? (a) : (-1.) )
  printf("# GYRATION %i %10.20g %10.20g %10.20g\n", i, m0, m1/NORM(m0), m2/NORM(m0));
}

return(0);
}



int gplist_corr(void)
{
int i, j, pk, pj, pk_dash;
int latest_j;
double corr_window_ave, corr, count;

#define POOL gplist_edge_pool_sorted
#define FETCH gplist_fetch_edge_sorted
#define LINK link_sorted
//#define POOL gplist_edge_pool_sequential
//#define FETCH gplist_fetch_edge_sequential

for (gplist_window_memory_pos=0; gplist_window_memory_pos<GPLIST_CORR_WINDOW_SIZE; gplist_window_memory_pos++) {
  MALLOC(gplist_window_memory[gplist_window_memory_pos], gplist_link_num);
  for (i=0; i<gplist_link_num; i++) {
    gplist_window_memory[gplist_window_memory_pos][i]=0.;
  }
}
gplist_window_memory_pos=0;
MALLOC(gplist_window_sum, gplist_link_num);
for (i=0; i<gplist_link_num; i++) 
  gplist_window_sum[i]=0.;


for (i=0; i<gplist_link_num; i++) {
  #warning "Maybe I should not be skipping no-neighbours."
  if (gplist_location[i].num_neighbours==0) continue;
  gplist_location[i].corr=0;
  /* Need to erase them all, because I don't know who I will encounter. */
  gplist_window_memory_pos++;
  /* XXX Need to initialise at leaset gplist_window_memory and gplist_window_sum! */
  for (j=0; j<gplist_link_num; j++) {
    gplist_window_memory[gplist_window_memory_pos%GPLIST_CORR_WINDOW_SIZE][j]=0;
    gplist_window_sum[j]-=gplist_window_memory[(gplist_window_memory_pos+1)%GPLIST_CORR_WINDOW_SIZE][j];
  }
  latest_j=0;
  for (j=0; j<gplist_link_num; j++) {
    if ((pj=FETCH(i, j))!=GPLIST_END_OF_LIST) {
      /* Supposedly, there is a count for i to j. We should find the same count from j to i. */
      if ((pk_dash=FETCH(j, i))==GPLIST_END_OF_LIST) { printf("# LINK_CLASH: %i %i (GPLIST_END_OF_LIST)\n", i, j); }
      else {
	if (POOL[pj].count!=POOL[pk_dash].count) { printf("# LINK_CLASH: %i %i (%g vs %g)\n", i, j, POOL[pj].count, POOL[pk_dash].count); }
	//printf("# COUNTS: %i %i: %g %g\n", i, j, POOL[pj].count, POOL[pk_dash].count);
      }
    }
    

    /* i is the index of the snippet we are looking at.
     * j is the index of one of its targets -- gplist_edge_pool_sequential[pj].count does not vanish.
     * pk runs through all links that i has
     * pk_dash is the corresponding link from j 
     */


    for (corr=0., pk=gplist_location[i].LINK; pk!=GPLIST_END_OF_LIST; pk=POOL[pk].next) {
      if ((pk_dash=FETCH(j, POOL[pk].target))==GPLIST_END_OF_LIST) continue;
      corr+= (POOL[pk].count * POOL[pk_dash].count );
    }
    if ((pj=FETCH(i, j))!=GPLIST_END_OF_LIST) {
      /* Supposedly, there is a count for i to j. We should find the same count from j to i. */
      if ((pk_dash=FETCH(j, i))==GPLIST_END_OF_LIST) { printf("# LINK_CLASH: %i %i (GPLIST_END_OF_LIST)\n", i, j); }
      else {
	if (POOL[pj].count!=POOL[pk_dash].count) { printf("# LINK_CLASH: %i %i (%g vs %g)\n", i, j, POOL[pj].count, POOL[pk_dash].count); }
	//printf("# COUNTS: %i %i: %g %g\n", i, j, POOL[pj].count, POOL[pk_dash].count);
      }
      POOL[pj].corr=corr;
      count=POOL[pj].count;
    } else count=0.;

    gplist_window_memory[gplist_window_memory_pos%GPLIST_CORR_WINDOW_SIZE][j]=corr;
    gplist_window_sum[j]+=corr;
    if (gplist_window_memory_pos>=GPLIST_CORR_WINDOW_SIZE-1) {
      /* In the sum are: GPLIST_CORR_WINDOW_SIZE-1, because I have just written one (the one that was deleted earlier) and I have subtracted one as well. 
       * At the very beginning, gplist_window_memory[0] will still be filled with 0s at the time when I land here for the first time. */
      corr_window_ave=gplist_window_sum[j]/((double)(GPLIST_CORR_WINDOW_SIZE-1));
    } else corr_window_ave=0.;
    latest_j=j;
    //printf("# gplist_window_sum=%g %i\n", gplist_window_sum[j], gplist_window_memory_pos);
    printf("# CORR: %i %i %g %g %g %g %g\n", i, j, corr, gplist_location[i].total_count, gplist_location[j].total_count, count, corr_window_ave);
  }
  while (latest_j<gplist_link_num-1) {
    latest_j++;
    if (gplist_window_memory_pos>=GPLIST_CORR_WINDOW_SIZE-1) {
      corr_window_ave=gplist_window_sum[latest_j]/((double)(GPLIST_CORR_WINDOW_SIZE-1));
    } else corr_window_ave=0.;
    printf("# CORR: %i %i %g %g %g %g %g\n", i, latest_j, 0., gplist_location[i].total_count, gplist_location[latest_j].total_count, 0., corr_window_ave);
  }
}
    
for (gplist_window_memory_pos=0; gplist_window_memory_pos<GPLIST_CORR_WINDOW_SIZE; gplist_window_memory_pos++) {
  free(gplist_window_memory[gplist_window_memory_pos]);
}
free(gplist_window_sum);
    
return(0);

}





int gplist_corr_count_weighted(void)
{
int i, j, pk, pj, pk_dash;
int latest_j;
double corr_window_ave;

#define POOL gplist_edge_pool_sorted
#define FETCH gplist_fetch_edge_sorted
#define LINK link_sorted
//#define POOL gplist_edge_pool_sequential
//#define FETCH gplist_fetch_edge_sequential

for (gplist_window_memory_pos=0; gplist_window_memory_pos<GPLIST_CORR_WINDOW_SIZE; gplist_window_memory_pos++) {
  MALLOC(gplist_window_memory[gplist_window_memory_pos], gplist_link_num);
  for (i=0; i<gplist_link_num; i++) {
    gplist_window_memory[gplist_window_memory_pos][i]=0.;
  }
}
gplist_window_memory_pos=0;
MALLOC(gplist_window_sum, gplist_link_num);
for (i=0; i<gplist_link_num; i++) 
  gplist_window_sum[i]=0.;


for (i=0; i<gplist_link_num; i++) {
  #warning "Maybe I should not be skipping no-neighbours."
  if (gplist_location[i].num_neighbours==0) continue;
  gplist_location[i].corr=0;
  /* Need to erase them all, because I don't know who I will encounter. */
  gplist_window_memory_pos++;
  /* XXX Need to initialise at leaset gplist_window_memory and gplist_window_sum! */
  for (j=0; j<gplist_link_num; j++) {
    gplist_window_memory[gplist_window_memory_pos%GPLIST_CORR_WINDOW_SIZE][j]=0;
    gplist_window_sum[j]-=gplist_window_memory[(gplist_window_memory_pos+1)%GPLIST_CORR_WINDOW_SIZE][j];
  }
  latest_j=0;
  for (pj=gplist_location[i].LINK; pj!=GPLIST_END_OF_LIST; pj=POOL[pj].next) {
    j=POOL[pj].target;
    while (latest_j<j-1) {
      latest_j++;
      if (gplist_window_memory_pos>=GPLIST_CORR_WINDOW_SIZE-1) {
        corr_window_ave=gplist_window_sum[latest_j]/((double)(GPLIST_CORR_WINDOW_SIZE-1));
      } else corr_window_ave=0.;
      printf("# CORR: %i %i %g %g %g %g %g\n", i, latest_j, 0., gplist_location[i].total_count, gplist_location[latest_j].total_count, 0., corr_window_ave);
    }


    {
    /* Supposedly, there is a count for i to j. We should find the same count from j to i. */
    if ((pk_dash=FETCH(j, i))==GPLIST_END_OF_LIST) { printf("# LINK_CLASH: %i %i (GPLIST_END_OF_LIST)\n", i, j); }
    else {
      if (POOL[pj].count!=POOL[pk_dash].count) { printf("# LINK_CLASH: %i %i (%g vs %g)\n", i, j, POOL[pj].count, POOL[pk_dash].count); }
      //printf("# COUNTS: %i %i: %g %g\n", i, j, POOL[pj].count, POOL[pk_dash].count);
    }
    }
    

    /* i is the index of the snippet we are looking at.
     * j is the index of one of its targets -- gplist_edge_pool_sequential[pj].count does not vanish.
     * pk runs through all links that i has
     * pk_dash is the corresponding link from j 
     */


    for (pk=gplist_location[i].LINK; pk!=GPLIST_END_OF_LIST; pk=POOL[pk].next) {
      if ((pk_dash=FETCH(j, POOL[pk].target))==GPLIST_END_OF_LIST) continue;
      POOL[pj].corr+= (POOL[pj].count*( POOL[pk].count * POOL[pk_dash].count ));
    }
    gplist_window_memory[gplist_window_memory_pos%GPLIST_CORR_WINDOW_SIZE][j]=POOL[pj].corr;
    gplist_window_sum[j]+=POOL[pj].corr;
    if (gplist_window_memory_pos>=GPLIST_CORR_WINDOW_SIZE-1) {
      /* In the sum are: GPLIST_CORR_WINDOW_SIZE-1, because I have just written one (the one that was deleted earlier) and I have subtracted one as well. 
       * At the very beginning, gplist_window_memory[0] will still be filled with 0s at the time when I land here for the first time. */
      corr_window_ave=gplist_window_sum[j]/((double)(GPLIST_CORR_WINDOW_SIZE-1));
    } else corr_window_ave=0.;
    latest_j=j;
    //printf("# gplist_window_sum=%g %i\n", gplist_window_sum[j], gplist_window_memory_pos);
    printf("# CORR: %i %i %g %g %g %g %g\n", i, j, POOL[pj].corr, gplist_location[i].total_count, gplist_location[j].total_count, POOL[pj].count, corr_window_ave);
  }
  while (latest_j<gplist_link_num-1) {
    latest_j++;
    if (gplist_window_memory_pos>=GPLIST_CORR_WINDOW_SIZE-1) {
      corr_window_ave=gplist_window_sum[latest_j]/((double)(GPLIST_CORR_WINDOW_SIZE-1));
    } else corr_window_ave=0.;
    printf("# CORR: %i %i %g %g %g %g %g\n", i, latest_j, 0., gplist_location[i].total_count, gplist_location[latest_j].total_count, 0., corr_window_ave);
  }
}
    
for (gplist_window_memory_pos=0; gplist_window_memory_pos<GPLIST_CORR_WINDOW_SIZE; gplist_window_memory_pos++) {
  free(gplist_window_memory[gplist_window_memory_pos]);
}
free(gplist_window_sum);
    
return(0);

}




int gplist_print_edges_sorted(void)
{
int source;
int p;

for (source=0; source<gplist_link_num; source++) 
  printf("#PRINT_META_DATA %5i: %i neighbours, pool starting at %i\n", source, gplist_location[source].num_neighbours, gplist_location[source].link_sorted);

for (source=0; source<gplist_link_num; source++) {
  printf("#PRINT_EDGES %5i:", source);
  for (p=gplist_location[source].link_sorted; (p!=GPLIST_END_OF_LIST); p=gplist_edge_pool_sorted[p].next) {
    printf("(T%i,C%g)", gplist_edge_pool_sorted[p].target, gplist_edge_pool_sorted[p].count);
  }
  printf("\n");
}
return(0);
}

int gplist_print_edges(void)
{
int source;
int p;


for (source=0; source<gplist_link_num; source++) {
  printf("#PRINT_EDGES %5i:", source);
  for (p=gplist_link[source]; (p!=GPLIST_END_OF_LIST); p=gplist_edge_pool_sequential[p].next) {
    printf("(T%i,C%g)", gplist_edge_pool_sequential[p].target, gplist_edge_pool_sequential[p].count);
  }
  printf("\n");
}
return(0);
}

int gplist_traversal(int visited) 
{
int ident=0;
int p, best_p=GPLIST_END_OF_LIST;

for (ident=0; ident<gplist_link_num; ident++) {
  gplist_visited[ident]=0;
}
  
for (ident=0; ;) {
  printf("#VISIT %i %i %i\n", visited, ident, (best_p==GPLIST_END_OF_LIST) ? 0 : 1);

  gplist_visited[ident]=visited;
  best_p=GPLIST_END_OF_LIST;
  for (p=gplist_link[ident]; (p!=GPLIST_END_OF_LIST); p=gplist_edge_pool_sequential[p].next) {
    /* find the strongest unused link */
    if (best_p==GPLIST_END_OF_LIST) {
      if (gplist_edge_pool_sequential[p].visited==0) 
        if (gplist_visited[gplist_edge_pool_sequential[p].target]==0) best_p=p;
    } else {
      if (gplist_edge_pool_sequential[p].visited==0) {
        if (gplist_edge_pool_sequential[best_p].count<gplist_edge_pool_sequential[p].count) 
	  if (gplist_visited[gplist_edge_pool_sequential[p].target]==0) best_p=p;
      }
    }
  }
  if (best_p==GPLIST_END_OF_LIST) {
    for (ident=0; ident<gplist_link_num; ident++) {
      if (gplist_visited[ident]==0) break;
    }
    if (ident==gplist_link_num) return(1);
    else continue;
  }
  ident=gplist_edge_pool_sequential[best_p].target;
  gplist_edge_pool_sequential[best_p].visited=visited;
}
return(0);
}
      

int gplist_write(char *filename)
{
FILE *file;

if ((file=fopen(filename, "w"))==NULL) {
  fprintf(stderr, "Failed to open file [%s] for writing, error %i::%s\n", filename, errno, strerror(errno));
  return(-1);
}

fwrite(&gplist_header, sizeof(gplist_header), 1, file);
fwrite(&gplist_link_num, sizeof(gplist_link_num), 1, file);
fwrite(gplist_link, sizeof(*gplist_link), gplist_link_num, file);
fwrite(gplist_visited, sizeof(*gplist_visited), gplist_link_num, file);
fwrite(&gplist_edge_pool_sequential_num, sizeof(gplist_edge_pool_sequential_num), 1, file);
fwrite(&gplist_edge_pool_sequential_count, sizeof(gplist_edge_pool_sequential_count), 1, file);
fwrite(gplist_edge_pool_sequential, sizeof(*gplist_edge_pool_sequential), gplist_edge_pool_sequential_count, file);

return(0);
}

int gplist_read(char *filename)
{
FILE *file;

if ((file=fopen(filename, "r"))==NULL) {
  fprintf(stderr, "Failed to open file [%s] for reading, error %i::%s\n", filename, errno, strerror(errno));
  return(-1);
}

fread(&gplist_header, sizeof(gplist_header), 1, file);
fread(&gplist_link_num, sizeof(gplist_link_num), 1, file);
MALLOC(gplist_link, gplist_link_num);
MALLOC(gplist_visited, gplist_link_num);
fread(gplist_link, sizeof(*gplist_link), gplist_link_num, file);
fread(gplist_visited, sizeof(*gplist_visited), gplist_link_num, file);

fread(&gplist_edge_pool_sequential_num, sizeof(gplist_edge_pool_sequential_num), 1, file);
fread(&gplist_edge_pool_sequential_count, sizeof(gplist_edge_pool_sequential_count), 1, file);
gplist_edge_pool_sequential_num=gplist_edge_pool_sequential_count;
MALLOC(gplist_edge_pool_sequential, gplist_edge_pool_sequential_count);
fread(gplist_edge_pool_sequential, sizeof(*gplist_edge_pool_sequential), gplist_edge_pool_sequential_count, file);

printf("# Info: Read %i snippets with %i edges.\n", gplist_link_num, gplist_edge_pool_sequential_count);

fclose(file);
return(0);
}



int gplist_prep_header(long int flags, char *command)
{
strcpy(gplist_header.magic, GPLIST_MAGIC);
strcpy(gplist_header.version, GPLIST_VERSION);
gplist_header.flags=flags;
gplist_header.flags |= ( (gplist_no_coarse_graining) ? (GPLIST_HDRFLG_NO_COARSE_GRAINING) : (GPLIST_HDRFLG_COARSE_GRAINING) );

gplist_header.link_num=gplist_link_num;
gplist_header.max_num_neighbours=gplist_max_num_neighbours;
if (flags & GPLIST_HDRFLG_SORTED) {
  gplist_header.edge_pool_count=gplist_edge_pool_sorted_count;
}
if (flags & GPLIST_HDRFLG_SEQUENTIAL) {
  gplist_header.edge_pool_count=gplist_edge_pool_sequential_count;
}
gplist_header.chr_w_length=gplist_chr_w_length;
gplist_header.Delta_w=gplist_Delta_w;
strncpy(gplist_header.generating_command, command, sizeof(gplist_header.generating_command));
gplist_header.generating_command[sizeof(gplist_header.generating_command)-1]=(char)0;
gplist_header.header_size=sizeof(gplist_header);
gplist_header.edge_size=sizeof(gplist_edge_strct);
gplist_header.location_size=sizeof(gplist_location_strct);
return(0);
}




int gplist_write_sorted(char *filename)
{
FILE *file;

if ((file=fopen(filename, "w"))==NULL) {
  fprintf(stderr, "Failed to open file [%s] for writing, error %i::%s\n", filename, errno, strerror(errno));
  return(-1);
}

fwrite(&gplist_header, sizeof(gplist_header), 1, file);
fwrite(gplist_location, sizeof(*gplist_location), gplist_link_num, file);
fwrite(gplist_edge_pool_sorted, sizeof(*gplist_edge_pool_sorted), gplist_edge_pool_sorted_count, file);

return(0);
}

int gplist_read_sorted(char *filename)
{
FILE *file;

printf("# Info: CVS Header of this file: $Header: /home/ma/p/pruess/.cvsroot/gplist/gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $\n");
printf("# Info: CVS/git ID of this file: $Id: gplist_util.c,v 1.2 2019/05/11 21:52:57 pruess Exp $\n");

if (filename[0]==(char)0) {
  file=stdin;
} else {
  if ((file=fopen(filename, "r"))==NULL) {
    fprintf(stderr, "Failed to open file [%s] for reading, error %i::%s\n", filename, errno, strerror(errno));
    return(-1);
  }
}
fread(&gplist_header, sizeof(gplist_header), 1, file);

if (strcmp(gplist_header.magic, GPLIST_MAGIC)) {
  fprintf(stderr, "GPLIST_MAGIC in header not found.\n");
  exit(EXIT_FAILURE);
}

gplist_link_num=gplist_header.link_num;
gplist_max_num_neighbours=gplist_header.max_num_neighbours;
gplist_edge_pool_sorted_count=gplist_header.edge_pool_count;

printf("# Info: Header of file [%s]\n", filename);
printf("# Info: Version: %s\n", gplist_header.version);
printf("# Info: link_num=%i\n", gplist_header.link_num);
printf("# Info: max_num_neighbours=%i\n", gplist_header.max_num_neighbours);
printf("# Info: edge_pool_count=%i\n", gplist_header.edge_pool_count);
printf("# Info: chr_w_length=%g\n", gplist_header.chr_w_length);
printf("# Info: Delta_w=%g\n", gplist_header.Delta_w);
printf("# Info: generating_command: %s\n", gplist_header.generating_command);
printf("# Info: header_size=%i versus %i\n", gplist_header.header_size, (int)sizeof(gplist_header));
printf("# Info: edge_size=%i versus %i\n", gplist_header.edge_size, (int)sizeof(gplist_edge_strct));
printf("# Info: location_size=%i versus %i\n", gplist_header.location_size, (int)sizeof(gplist_location_strct));

if (gplist_header.link_num!= (((int)(gplist_header.chr_w_length/gplist_header.Delta_w))+1) ) {
  fprintf(stderr, "gplist_header.link_num=%i != ((int)gplist_header.chr_w_length/gplist_header.Delta_w)+1)=%i\n", gplist_header.link_num, ((int)(gplist_header.chr_w_length/gplist_header.Delta_w))+1);
  exit(EXIT_FAILURE);
}

MALLOC(gplist_location, gplist_link_num);
MALLOC(gplist_edge_pool_sorted, gplist_edge_pool_sorted_count);

fread(gplist_location, sizeof(*gplist_location), gplist_link_num, file);
fread(gplist_edge_pool_sorted, sizeof(*gplist_edge_pool_sorted), gplist_edge_pool_sorted_count, file);

printf("# Info: Read %i snippets with %i edges.\n", gplist_link_num, gplist_edge_pool_sorted_count);
if (file!=stdin) fclose(file);

return(0);
}















