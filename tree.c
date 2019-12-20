/*  File: tree.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 19 10:53 2019 (rd109)
 * Created: Sat Nov 16 10:07:06 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "tree.h"
#include <math.h>

void  treeNodeDestroy (TreeNode *n) // recursive destroy
{
  if (n)
    { treeNodeDestroy (n->left) ;
      treeNodeDestroy (n->right) ;
      if (n->isRoot && n->parent) treeNodeDestroy (n->parent) ;
      if (n->name) free (n->name) ;
      free (n) ;
    }
}

static int addElement (Tree *t, TreeNode *n, int kParent)
{
  static char nameBuf[32] ;
  assert ((n->left && n->right) || (!n->left && !n->right)) ; // check well-formed
  int kSelf = arrayMax(t->a) ;
  { char *name ; int i ;
    if (n->name) name = n->name ; else { sprintf (nameBuf, "_%d", kSelf) ; name = nameBuf ; }
    if (!dictAdd (t->nameDict, name, &i)) die ("duplicate tree node name %s", name) ;
    assert (i == kSelf) ;
  }
  TreeElement *e = arrayp(t->a, arrayMax(t->a), TreeElement) ;
  e->length = n->length ;
  if (e->length < 0)
    { fprintf (stderr, "negative length %.4g fixed for node %d\n", e->length, kSelf) ;
      e->length = -e->length ;
    }

  // next the recursive section - need to use local eLeft, eRight etc. since e is volatile
  int eLeft = 0, eRight = 0, eParent = 0 ;
  if (n->left) eLeft = addElement (t, n->left, kSelf) ;
  if (n->right) eRight = addElement (t, n->right, kSelf) ;
  if (n->isRoot && n->parent) eParent = addElement (t, n->parent, kSelf) ;
  e = arrayp(t->a, kSelf, TreeElement) ; // reassign e, in case memory for t->a moved
  e->left = eLeft ;
  e->right = eRight ;
  e->parent = eParent ? eParent : kParent ;

  return kSelf ;
}

Tree *treeCreate (TreeNode *n)
{
  Tree *t = new0 (1, Tree) ;
  t->a = arrayCreate (1024, TreeElement) ;
  t->nameDict = dictCreate (1024) ;
  addElement (t, n, 0) ;

  return t ;
}

void treeDestroy (Tree *t)
{
  arrayDestroy (t->a) ;
  dictDestroy (t->nameDict) ;
}

/*****************************************/

static void reportHeightStats (Tree *t)
{
  int i, n = 0 ;
  double *rootDist = new0 (arrayMax(t->a), double) ;
  double sum = 0., sum2 = 0. ;
  for (i = 0 ; i < arrayMax(t->a) ; ++i) // pre-order
    { TreeElement *e = arrp(t->a, i, TreeElement) ;
      if (e->parent >= 0) rootDist[i] = rootDist[e->parent] + e->length ;
      if (!e->left) { ++n ; sum += rootDist[i] ; sum2 += rootDist[i]*rootDist[i] ; } // a leaf
    }
  printf ("n, mean, sd root distance of leaves: %d %.4g %.4g\n",
	  n, sum/n, sqrt((sum2-(sum*sum)/n)/(n-1))) ;
}

void treeBalance (Tree *t)
{
  reportHeightStats (t) ;
  int i ;
  double *height = new0 (arrayMax(t->a), double) ; // height at top of edge
  double *factor = new0 (arrayMax(t->a), double) ;
  for (i = arrayMax(t->a) ; i-- ; ) // first build up in post-order
    { TreeElement *e = arrp(t->a, i, TreeElement) ;
      height[i] = e->length ;
      if (e->left)
	{ double mean = (height[e->left] + height[e->right]) / 2. ;
	  height[i] += mean ;
	  factor[e->left] = mean / height[e->left] ;
	  factor[e->right] = mean / height[e->right] ;
	}
    }
  for (i = 0 ; i < arrayMax(t->a) ; ++i) // now go back down and fix in pre-order
    if (factor[i])
      { TreeElement *e = arrp(t->a, i, TreeElement) ;
	e->length *= factor[i] ;
	if (e->left)
	  { factor[e->left] *= factor[i] ;
	    factor[e->right] *= factor[i] ;
	  }
      }
  free (height) ; free (factor) ;
  reportHeightStats (t) ;
}

/*****************************************/

typedef struct { double pFlip, pStick ; } TreeEdge ;

Array treeBuildEdges (Tree *t, double rate)
{
  Array a = arrayCreate (arrayMax(t->a), TreeEdge) ;
  int i ;
  for (i = arrayMax(t->a) ; i-- ; )
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      TreeEdge *d = arrayp(a, i, TreeEdge) ;
      if (e->length)
	{ d->pFlip = (1. - exp(-rate*e->length)) / 2. ;
	  d->pStick = (1. + exp(-rate*e->length)) / 2. ;
	}
    }
  return a ;
}

static inline void addScore (LogLikelihood *a, LogLikelihood *b, TreeEdge *e)
{
  if (b->s0 == b->s1) // e.g. at empty leaves
    { a->s0 += b->s0 ; a->s1 += b->s1 ; }
  else 
    { double ratio = exp (b->s1 - b->s0) ;
      a->s0 += b->s0 + log (e->pStick + e->pFlip * ratio) ;
      a->s1 += b->s1 + log (e->pStick + e->pFlip / ratio) ;
    }
  assert (!isnan(a->s0) && !isnan(a->s1)) ;
}

TreeScore *treeBuildScores (Tree *t, char *gt, int* tree2vcf, Array edges)
{
  static float LEAF_MATCH = -0.00010005, LEAF_MISMATCH = -9.21034037 ; // log(0.9999), log(0.0001)
  int i ;
  TreeScore *scores = new0 (arrayMax(t->a), TreeScore) ;
  // first make scores below by running from bottom up the tree
  for (i = arrayMax(t->a) ; i-- ; ) // NB important i is descending for post-order
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      LogLikelihood *si = &(scores[i].below) ;
      if (e->left) // an internal node
	{ addScore (si, &scores[e->left].below, arrp(edges,e->left,TreeEdge)) ;
	  addScore (si, &scores[e->right].below, arrp(edges,e->right,TreeEdge)) ;
	  if (!i && e->parent) // third branch in unrooted tree
	    addScore (si, &scores[e->parent].below, arrp(edges,e->parent,TreeEdge)) ;
	}
      else // a leaf
	switch (gt[tree2vcf[i]])
	  {
	  case 3: si->s0 = LEAF_MISMATCH ; si->s1 = LEAF_MATCH ; break ;
	  case 2: si->s0 = LEAF_MATCH ; si->s1 = LEAF_MISMATCH ; break ;
	  case 1: si->s0 = si->s1 = 0. ; break ;
	  default: // die ("unknown gt %d at position %d", gt[tree2vcf[i]], i) ;
	    fprintf (stderr, "bad gt %d at position %d %d\n", gt[tree2vcf[i]], i, tree2vcf[i]) ;
	  }
    }
  // then above, running from top down the tree
  for (i = 0 ; i < arrayMax(t->a) ; i++ ) // now i is ascending order for pre-order
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      if (!e->left) continue ; // make for the children of i
      LogLikelihood *sl = &(scores[e->left].above), *sr = &(scores[e->right].above) ;
      addScore (sl, &(scores[e->right].below), arrp(edges,e->right,TreeEdge)) ; // NB below
      addScore (sr, &(scores[e->left].below), arrp(edges,e->left,TreeEdge)) ;
      if (i) // not the root
	{ addScore (sl, &(scores[i].above), arrp(edges,i,TreeEdge)) ;
	  addScore (sr, &(scores[i].above), arrp(edges,i,TreeEdge)) ;
	}
      else if (e->parent) // third branch in unrooted tree 
	{ addScore (sl, &(scores[e->parent].below), arrp(edges,e->parent,TreeEdge)) ;
	  addScore (sr, &(scores[e->parent].below), arrp(edges,e->parent,TreeEdge)) ;
	  LogLikelihood *sp = &(scores[e->parent].above) ;
	  addScore (sp, &(scores[e->right].below), arrp(edges,e->right,TreeEdge)) ; // NB below
	  addScore (sp, &(scores[e->left].below), arrp(edges,e->left,TreeEdge)) ;
	}
    }
  // as currently coded this is a bit inefficient - I add the up-scores twice, once when
  // making below, and again when making above.  I could cache them with a bit more complexity.
  return scores ;
}

/***********************************************************************/

