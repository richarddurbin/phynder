/*  File: tree.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 30 12:29 2019 (rd109)
 * Created: Fri Nov 15 23:58:52 2019 (rd109)
 *-------------------------------------------------------------------
 */

#ifndef TREE_DEF
#define TREE_DEF

#include "utils.h"

typedef struct TreeNodeStruct {
  struct TreeNodeStruct *parent, *left, *right ;
  char *name ;
  double length ;
  BOOL isRoot ;
  BOOL mark ;
} TreeNode ;

typedef struct {
  int parent, left, right ;	/* indices in array */
  double length ;
} TreeElement ;

typedef struct {
  Array a ; // pre-order array of TreeElement, access in reverse for post-order
  DICT  *nameDict ;
} Tree ;

typedef struct { float s0, s1 ; } LogLikelihood ; // 's' for "score"

typedef struct
{ LogLikelihood below ; // log likelihoods of data below node given 0 or 1 at bottom of edge
  LogLikelihood above ; // log likelihoods of data above node given 0 or 1 at top of edge
} TreeScore ;

void  treeNodeDestroy (TreeNode *n) ; // recursive destroy
Tree *treeCreate (TreeNode *n) ;
void  treeDestroy (Tree *t) ;
Array treeRateBuild (Tree *t, double rate) ;
void  treeBalance (Tree *t) ;
Array treeBuildEdges (Tree *t, double rate) ;
TreeScore *treeBuildScores (Tree *t, char *gt, int* tree2vcf, Array edges) ;

#endif

/****************/
