/*  File: tree.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 10 13:59 2020 (rd109)
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

void  treeNodeDestroy (TreeNode *n) ; // recursive destroy
Tree *treeCreate (TreeNode *n) ;
void  treeDestroy (Tree *t) ;
Array treeRateBuild (Tree *t, double rate) ;
void  treeBalance (Tree *t) ;
Array treeBuildEdges (Tree *t, double rate, double *worst) ;
LogLikelihood *treeBuildScores (Tree *t, char *gt, int* tree2vcf, Array edges, int calcMode) ;

#endif

/****************/
