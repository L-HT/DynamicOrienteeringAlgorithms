/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef  __LINKERN_H
#define  __LINKERN_H

#include "LINKERN_support_data.h"
#include "LINKERN_support_flip_two.h"
#include "LINKERN_support_CCrandstate_cc.h"
#include "HelperFunctions.h"

#define CC_LK_RANDOM_KICK    (0)
#define CC_LK_GEOMETRIC_KICK (1)
#define CC_LK_CLOSE_KICK     (2)
#define CC_LK_WALK_KICK      (3)

#define MAK_MORTON
#undef  FULL_MAK_MORTON
#undef  NODE_INSERTIONS

#ifdef CC_ATTRIBUTE
#define CC_UNUSED __attribute__ ((unused))
#else
#define CC_UNUSED
#endif


////////////////////////////////////////////////


/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#define CC_SAFE_MALLOC(nnum,type)                                          \
    (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                             \
    CCutil_freerus ((void *) (object));                                    \
    object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
    if ((object)) CC_FREE ((object),type);                                 \
}

#define CC_PTRWORLD_ALLOC_ROUTINE(type, ptr_alloc_r, ptr_bulkalloc_r)        \
                                                                             \
static int ptr_bulkalloc_r (CCptrworld *world, int nalloc)                   \
{                                                                            \
    CCbigchunkptr *bp;                                                       \
    int i;                                                                   \
    int count = CC_BIGCHUNK / sizeof ( type );                               \
    type *p;                                                                 \
                                                                             \
    while (nalloc > 0) {                                                     \
        bp = CCutil_bigchunkalloc ();                                        \
        if (bp == (CCbigchunkptr *) NULL) {                                  \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return 1;                                                        \
        }                                                                    \
        bp->next = world->chunklist ;                                        \
        world->chunklist = bp;                                               \
                                                                             \
        p = ( type * ) bp->this_one;                                         \
        for (i=count-2; i>=0; i--) {                                         \
            p[i].next = &p[i+1];                                             \
        }                                                                    \
        p[count - 1].next = (type *) world->freelist;                        \
        world->freelist = (void *) p;                                        \
        nalloc -= count;                                                     \
    }                                                                        \
    return 0;                                                                \
}                                                                            \
                                                                             \
static type *ptr_alloc_r (CCptrworld *world)                                 \
{                                                                            \
    type *p;                                                                 \
                                                                             \
    if (world->freelist == (void *) NULL) {                                  \
        if (ptr_bulkalloc_r (world, 1)) {                                    \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return ( type * ) NULL;                                          \
        }                                                                    \
    }                                                                        \
    p = (type *) world->freelist ;                                           \
    world->freelist = (void *) p->next;                                      \
                                                                             \
    return p;                                                                \
}

#define CC_PTRWORLD_FREE_ROUTINE(type, ptr_free_r)                           \
                                                                             \
static void ptr_free_r (CCptrworld *world, type *p)                          \
{                                                                            \
    p->next = (type *) world->freelist ;                                     \
    world->freelist = (void *) p;                                            \
}

#define CC_PTRWORLD_LISTADD_ROUTINE(type, entrytype, ptr_listadd_r, ptr_alloc_r) \
                                                                             \
static int ptr_listadd_r (type **list, entrytype x, CCptrworld *world)       \
{                                                                            \
    if (list != (type **) NULL) {                                            \
        type *p = ptr_alloc_r (world);                                       \
                                                                             \
        if (p == (type *) NULL) {                                            \
            fprintf (stderr, "ptr list add failed\n");                       \
            return 1;                                                        \
        }                                                                    \
        p->this = x;                                                         \
        p->next = *list;                                                     \
        *list = p;                                                           \
    }                                                                        \
    return 0;                                                                \
}

#define CC_PTRWORLD_LISTFREE_ROUTINE(type, ptr_listfree_r, ptr_free_r)       \
                                                                             \
static void ptr_listfree_r (CCptrworld *world, type *p)                      \
{                                                                            \
    type *next;                                                              \
                                                                             \
    while (p != (type *) NULL) {                                             \
        next = p->next;                                                      \
        ptr_free_r (world, p);                                               \
        p = next;                                                            \
    }                                                                        \
}

#define CC_PTRWORLD_LEAKS_ROUTINE(type, ptr_leaks_r, field, fieldtype)       \
                                                                             \
static int ptr_leaks_r (CCptrworld *world, int *total, int *onlist)          \
{                                                                            \
    int count = CC_BIGCHUNK / sizeof ( type );                               \
    int duplicates = 0;                                                      \
    type * p;                                                                \
    CCbigchunkptr *bp;                                                       \
                                                                             \
    *total = 0;                                                              \
    *onlist = 0;                                                             \
                                                                             \
    for (bp = world->chunklist ; bp; bp = bp->next)                          \
        (*total) += count;                                                   \
                                                                             \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        (*onlist)++;                                                         \
        p-> field = ( fieldtype ) 0;                                         \
    }                                                                        \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        if ((unsigned long) p-> field == (unsigned long) (size_t) 1)                           \
            duplicates++;                                                    \
        else                                                                 \
            p-> field = ( fieldtype ) (size_t) 1;                            \
    }                                                                        \
    if (duplicates) {                                                        \
        fprintf (stderr, "WARNING: %d duplicates on ptr free list \n",       \
                 duplicates);                                                \
    }                                                                        \
    return *total - *onlist;                                                 \
}

#define CC_PTRWORLD_ROUTINES(type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r) \
CC_PTRWORLD_ALLOC_ROUTINE (type, ptr_alloc_r, ptr_bulkalloc_r)               \
CC_PTRWORLD_FREE_ROUTINE (type, ptr_free_r)

#define CC_PTRWORLD_LIST_ROUTINES(type, entrytype, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r, ptr_listadd_r, ptr_listfree_r) \
CC_PTRWORLD_ROUTINES (type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r)        \
CC_PTRWORLD_LISTADD_ROUTINE (type, entrytype, ptr_listadd_r, ptr_alloc_r)    \
CC_PTRWORLD_LISTFREE_ROUTINE (type, ptr_listfree_r, ptr_free_r)

#define CC_BIGCHUNK ((int) ((1<<16) - sizeof (CCbigchunkptr) - 16))

struct CCbigchunk;

typedef struct CCbigchunkptr {
    void                 *this_one;
    struct CCbigchunk    *this_chunk;
    struct CCbigchunkptr *next;
} CCbigchunkptr;


typedef struct CCptrworld {
    int refcount;
    void *freelist;
    CCbigchunkptr *chunklist;
} CCptrworld;


#define CC_BIGCHUNK ((int) ((1<<16) - sizeof (CCbigchunkptr) - 16))


typedef struct CCbigchunk {
    char space[CC_BIGCHUNK];
    CCbigchunkptr ptr;
} CCbigchunk;

typedef struct CCdheap {
    double  *key;
    int     *entry;
    int     *loc;
    int     total_space;
    int     size;
} CCdheap;

void CCptrworld_init (CCptrworld *world);

//////////////////////////////////////////////////

typedef struct edge {
    int other;
    int weight;
} edge;

typedef struct edgelook {
    struct edgelook *next;
    int other;
    int diff;
    int over;
    int seq;
    int side;
#ifdef MAK_MORTON
    int mm;
#endif
#ifdef NODE_INSERTIONS
    int ni;
    int under;
#endif
} edgelook;

typedef struct intptr {
    int this_2;
    struct intptr *next;
} intptr;

typedef struct flippair {
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack {
    flippair *stack;
    int counter;
    int max;
} flipstack;

typedef struct graph {
    edge **goodlist;
    edge *edgespace;
    int  *degree;
    int  *weirdmark;
    int   weirdmagic;
    int   ncount;
    CCrandstate *rstate;
} graph;

typedef struct distobj {
    compass_data *dat;
    int       *cacheval;
    int       *cacheind;
    int        cacheM;
} distobj;

typedef struct adddel {
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct aqueue {
    char *active;
    intptr *active_queue;
    intptr *bottom_active_queue;
    CCdheap *h;
} aqueue;

typedef struct CCkdtree{

} CCkdtree;


int
    CClinkern_tour (int ncount, compass_data *dat, int ecount,
        int *elist, int stallcount, int repeatcount, int *incycle,
        int *outcycle, double *val, int silent, double time_bound,
        double length_bound, char *saveit_name, int kicktype,
        CCrandstate *rstate, AdditionalLogData* logData),
    CClinkern_path (int ncount, compass_data *dat, int ecount,
        int *elist, int nkicks, int *inpath, int *outpath, double *val,
        int silent, CCrandstate *rstate),
    CClinkern_fixed (int ncount, compass_data *dat, int ecount, int *elist,
        int nkicks, int *incycle, int *outcycle, double *val, int fcount,
        int *flist, int silent, CCrandstate *rstate);

CCbigchunkptr *CCutil_bigchunkalloc (void);
void CCutil_bigchunkfree (CCbigchunkptr *bp);
void CCptrworld_init (CCptrworld *world);
void CCutil_dheap_free (CCdheap *h);
void CCptrworld_delete (CCptrworld *world);
double CCutil_zeit (void);
int CCutil_dat_edgelen (int i, int j, compass_data *dat);
/////
//CCbigchunkptr *CCutil_bigchunkalloc (void)
//{
    //CCbigchunk *p = CC_SAFE_MALLOC (1, CCbigchunk);

    //if (p == (CCbigchunk *) NULL) {
        //fprintf (stderr, "Out of memory in CCutil_bigchunkalloc\n");
        //return (CCbigchunkptr *) NULL;
    //}
    //p->ptr.this_chunk = p;
    //p->ptr.this_one = (void *) p->space;
    //return &(p->ptr);
//}



//void CCutil_bigchunkfree (CCbigchunkptr *bp)
//{
    ///* This copy is necessary since CC_FREE zeros its first argument */
    //CCbigchunk *p = bp->this_chunk;

    //CC_FREE (p, CCbigchunk);
//}

//void CCptrworld_init (CCptrworld *world)
//{
    //world->refcount = 1;
    //world->freelist = (void *) NULL;
    //world->chunklist = (CCbigchunkptr *) NULL;
//}

//void CCutil_dheap_free (CCdheap *h)
//{
    //CC_IFFREE (h->entry, int);
    //CC_IFFREE (h->loc, int);
    //CC_IFFREE (h->key, double);
//}

//void CCptrworld_delete (CCptrworld *world)
//{
    //world->refcount--;
    //if (world->refcount <= 0) {
        //CCbigchunkptr *bp, *bpnext;

        //for (bp = world->chunklist ; bp; bp = bpnext) {
            //bpnext = bp->next;
            //CCutil_bigchunkfree (bp);
        //}
        //world->chunklist = (CCbigchunkptr *) NULL;
        //world->freelist = (void *) NULL;
        //world->refcount = 0;
    //}
//}
//double CCutil_zeit (void)
//{
    //return 0.0;
//}
//int CCutil_dat_edgelen (int i, int j, compass_data *dat)
//{
    //if (dat->ndepot) {
        //if (i >= dat->orig_ncount) {
            //return dat->depotcost[j];
        //} else if (j >= dat->orig_ncount) {
            //return dat->depotcost[i];
        //}
    //}
    //return (dat->edgelen)(i, j, dat);
//}


#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CC_OURABS(a) (((a) >= 0) ? (a) : -(a))


#endif  /* __LINKERN_H */



/****************************************************************************/
/*                                                                          */
/*                  FLIPPER HEADER (TWO-LIST)                               */
/*                                                                          */
/****************************************************************************/


#ifndef __FLIPPER_H
#define __FLIPPER_H

typedef struct CClk_parentnode {
    struct CClk_parentnode *adj[2];
    struct CClk_childnode  *ends[2];
    int                     size;
    int                     id;
    int                     rev;
} CClk_parentnode;

typedef struct CClk_childnode {
    struct CClk_parentnode *parent;
    struct CClk_childnode  *adj[2];
    int                     id;
    int                     name;
} CClk_childnode;

typedef struct CClk_flipper {
    CClk_parentnode        *parents;
    CClk_childnode         *children;
    int                     reversed;
    int                     nsegments;
    int                     groupsize;
    int                     split_cutoff;
} CClk_flipper;



int
    CClinkern_flipper_init (CClk_flipper *f, int ncount, int *cyc),
    CClinkern_flipper_next (CClk_flipper *f, int x),
    CClinkern_flipper_prev (CClk_flipper *f, int x),
    CClinkern_flipper_sequence (CClk_flipper *f, int x, int y, int z);
void
    CClinkern_flipper_flip (CClk_flipper *F, int x, int y),
    CClinkern_flipper_cycle (CClk_flipper *F, int *x),
    CClinkern_flipper_finish (CClk_flipper *F);

#endif  /* __FLIPPER_H */

#ifndef ADDITIONAL_SUPPORT_FUNCTIONS
#define ADDITIONAL_SUPPORT_FUNCTIONS

int random_four_swap (graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
       CCkdtree *kdt, int *delta, int kicktype, flipstack *win,
       flipstack *fstack, CCptrworld *intptr_world, CCrandstate *rstate);

int find_geometric_four (graph *G, distobj *D, CClk_flipper *F,
        CCkdtree *kdt, int *t1, int *t2, int *t3, int *t4, int *t5, int *t6,
        int *t7, int *t8, CCrandstate *rstate);

#endif
