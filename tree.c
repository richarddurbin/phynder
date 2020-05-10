/*  File: tree.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 10 13:58 2020 (rd109)
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

Array treeBuildEdges (Tree *t, double rate, double *worst)
{
  Array a = arrayCreate (arrayMax(t->a), TreeEdge) ;
  int i ;
  double shortest = 0. ;
  for (i = arrayMax(t->a) ; i-- ; )
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      TreeEdge *d = arrayp(a, i, TreeEdge) ;
      if (e->length)
	{ d->pFlip = (1. - exp(-rate*e->length)) / 2. ;
	  d->pStick = (1. + exp(-rate*e->length)) / 2. ;
	  if (!shortest || e->length < shortest)
	    { *worst = log(d->pFlip) ;
	      shortest = e->length ;
	    }
	}
    }
  return a ;
}

static inline void addEdge (LogLikelihood *a, LogLikelihood *b, TreeEdge *e)
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

static LogLikelihood *treeBuildBelow (Tree *t, char *gt, int* tree2vcf, Array edges)
{
  // below[i].s0 = p(data below edge i | 0 at bottom of edge i)
  static float LEAF_MATCH = -0.00010005, LEAF_MISMATCH = -9.21034037 ; // log(0.9999), log(0.0001)
  LogLikelihood *below = new0 (arrayMax(t->a), LogLikelihood) ;
  int i ;
  for (i = arrayMax(t->a) ; i-- ; ) // NB important i is descending for post-order
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      LogLikelihood *bi = &below[i] ;
      if (e->left) // an internal node
	{ addEdge (bi, &below[e->left], arrp(edges,e->left,TreeEdge)) ;
	  addEdge (bi, &below[e->right], arrp(edges,e->right,TreeEdge)) ;
	  if (!i && e->parent) // third branch in unrooted tree
	    addEdge (bi, &below[e->parent], arrp(edges,e->parent,TreeEdge)) ;
	}
      else // a leaf
	switch (gt[tree2vcf[i]])
	  {
	  case 3: bi->s0 = LEAF_MISMATCH ; bi->s1 = LEAF_MATCH ; break ;
	  case 2: bi->s0 = LEAF_MATCH ; bi->s1 = LEAF_MISMATCH ; break ;
	  case 1: bi->s0 = bi->s1 = 0. ; break ;
	  default: // die ("unknown gt %d at position %d", gt[tree2vcf[i]], i) ;
	    fprintf (stderr, "bad gt %d at position %d %d\n", gt[tree2vcf[i]], i, tree2vcf[i]) ;
	  }
    }
  return below ;
}

static LogLikelihood *treeBuildAbove (Tree *t, LogLikelihood *below, Array edges)
{
  // above[i].s0 = p(data not below edge i | 0 at top of edge i)
  LogLikelihood *above = new0 (arrayMax(t->a), LogLikelihood) ;
  int i ;
  for (i = 0 ; i < arrayMax(t->a) ; i++ ) // now i is ascending order for pre-order
    { TreeElement *e = arrp (t->a, i, TreeElement) ;
      if (!e->left) continue ; // make for the children of i
      LogLikelihood *al = &above[e->left], *ar = &above[e->right] ;
      addEdge (al, &below[e->right], arrp(edges,e->right,TreeEdge)) ; // NB below
      addEdge (ar, &below[e->left], arrp(edges,e->left,TreeEdge)) ;
      if (i) // not the root
	{ addEdge (al, &above[i], arrp(edges,i,TreeEdge)) ;
	  addEdge (ar, &above[i], arrp(edges,i,TreeEdge)) ;
	}
      else if (e->parent) // third branch in unrooted tree 
	{ addEdge (al, &below[e->parent], arrp(edges,e->parent,TreeEdge)) ;
	  addEdge (ar, &below[e->parent], arrp(edges,e->parent,TreeEdge)) ;
	  LogLikelihood *ap = &above[e->parent] ;
	  addEdge (ap, &below[e->right], arrp(edges,e->right,TreeEdge)) ; // NB below
	  addEdge (ap, &below[e->left], arrp(edges,e->left,TreeEdge)) ;
	}
    }
  return above ;
}

LogLikelihood *treeBuildScores (Tree *t, char *gt, int* tree2vcf, Array edges, int calcMode)
{
  LogLikelihood *below = treeBuildBelow (t, gt, tree2vcf, edges) ;
  LogLikelihood *above = treeBuildAbove (t, below, edges) ;
  // as currently coded this is a bit inefficient - I add the up-scores twice, once when
  // making below, and again when making above.  I could cache them with a bit more complexity.
  LogLikelihood *scores = new0 (arrayMax(t->a), LogLikelihood) ;

  // scores[0] is the LL at the root
  scores[0].s0 = below[0].s0 ; scores[0].s1 = below[0].s1 ;
  int i ;
  switch (calcMode)
    {
    case -1: // LL switch on this edge
      for (i = 1 ; i < arrayMax(t->a) ; ++i)
	{ scores[i].s0 = above[i].s0 + below[i].s1 ;
	  scores[i].s1 = above[i].s1 + below[i].s0 ;
	}
      break ;
    case 0: // LL unchanged match
      for (i = 1 ; i < arrayMax(t->a) ; ++i)
	{ scores[i].s0 = above[i].s0 + below[i].s0 ;
	  scores[i].s1 = above[i].s1 + below[i].s1 ;
	}
      break ;
    case 1: // LL not unchanged mismatch - not written yet, so for now as case 0
      for (i = 1 ; i < arrayMax(t->a) ; ++i)
	{ scores[i].s0 = above[i].s0 + below[i].s0 ;
	  scores[i].s1 = above[i].s1 + below[i].s1 ;
	}
      break ;
    default: die ("unknown calc mode %d", calcMode) ;
    }

  free (below) ;
  free (above) ;
  return scores ;
}

/***********************************************************************/
