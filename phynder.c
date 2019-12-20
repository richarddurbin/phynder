/*  File: phynder.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 19 18:35 2019 (rd109)
 * Created: Sun Nov 17 19:52:20 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "tree.h"
#include "newick.h"
#include "vcf.h"
#include <ctype.h>
#include <math.h>

static double transitionRate = 1.33 ; // expect ratio 4:1 for ts:tv rates
static double transversionRate = 0.33 ; // these values give mean 1.0 for 2:1 ts:tv events

static BOOL isVerbose = FALSE ;

static int baseMap[128] ;
#define POS_HASH(s) (((s)->pos << 4) + (baseMap[(s)->ref] << 2) + baseMap[(s)->alt])

void usage (void)
{
  fprintf (stderr, "usage: phynder [arguments/options]\n") ;
  fprintf (stderr, "  -t <newick tree>\n") ;
  fprintf (stderr, "  -v <vcf for tree>        must follow -t\n") ;
  fprintf (stderr, "  -q <query vcf>           must follow -t and -v\n") ;
  fprintf (stderr, "  -ts <transition rate>    [%.4f]\n", transitionRate) ;
  fprintf (stderr, "  -tv <transversion rate>  [%.4f]\n", transversionRate) ;
  fprintf (stderr, "  -B <branch file name>    output branch positions of tree variants\n") ;
  fprintf (stderr, "  -h                       print this message\n") ;
  fprintf (stderr, "  -V                       verbose - print extra info\n") ;
  exit (0) ;
}

int main (int argc, char *argv[])
{
  FILE *f ;
  Tree *t = 0 ;
  Vcf *vt = 0, *vq = 0 ;
  FILE *fB = 0 ;
 
  timeUpdate (0) ;

  baseMap['a'] = baseMap['A'] = 0 ;
  baseMap['c'] = baseMap['C'] = 1 ;
  baseMap['g'] = baseMap['G'] = 2 ;
  baseMap['t'] = baseMap['T'] = 3 ;

  --argc ; ++argv ; // absorb executable name
  if (!argc || !strcmp (*argv, "-h")) usage () ;
  while (argc)
    if (!strcmp (*argv, "-t") && argc > 1)
      { if (!(f = fopen (argv[1], "r"))) die ("failed to open newick tree file %s", argv[1]) ;
	TreeNode *n = readBinaryNewickTree (f) ;
	fclose (f) ;
	t = treeCreate (n) ;
	printf ("read tree with %d nodes and %d leaves from %s\n",
		arrayMax(t->a), (arrayMax(t->a)+1)/2, argv[1]) ;
	treeNodeDestroy (n) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-v"))
      { vt = vcfRead (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-q"))
      { vq = vcfRead (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-ts") && argc > 1)
      { transitionRate = atof (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-tv") && argc > 1)
      { transversionRate = atof (argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-B") && argc > 1)
      { if (!(fB = fopen (argv[1], "w"))) die ("failed to write branch file %s", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-V"))
      { isVerbose = TRUE ;
	--argc ; ++argv ;
      }
    else
      die ("unrecognized command line - run with no arguments for usage") ;

  if (!t) die ("must specify a tree with -t argument") ;

  treeBalance (t) ;
  timeUpdate (stdout) ;
  
  Array transitionEdges = treeBuildEdges (t, transitionRate) ;
  Array transversionEdges = treeBuildEdges (t, transversionRate) ;

  // build map from tree vcf to tree leaf nodes
  if (!vt) die ("must give a vcf for the tree") ;
  int *tree2vcf = new0(arrayMax(t->a), int) ;
  BOOL *leafFound = new0(arrayMax(t->a), BOOL) ;
  int i, j, k ;
  for (j = 0 ; j < dictMax(vt->samples) ; ++j)
    if (dictFind (t->nameDict, dictName (vt->samples,j), &i))
      { if (arrp(t->a,i,TreeElement)->left) continue ; // only interested in leaves
	tree2vcf[i] = j ;
	leafFound[i] = TRUE ;
      }
  for (i = 0 ; i < arrayMax(t->a) ; ++i)
    if (!arrp(t->a,i,TreeElement)->left && !leafFound[i])
      { char *name = dictName (t->nameDict, i) ;
	if (!dictFind (vt->samples, name, &j)) die ("missing sample %s in tree vcf\n", name) ;
	else die ("problem connecting tree leaf %s i %d to vcf sample j %d", name, i, j) ;
      }
  free (leafFound) ;
  printf ("built index from tree to vcf\n") ;
  timeUpdate (stdout) ;
  
  // build inside and outside log likeliHoods for each site for each node in tree
  // indirection via siteIndex so that info for sites with identical genotypes is shared
  TreeScore **scores = new0(arrayMax(vt->sites), TreeScore*) ; // more than needed but cheap
  int *siteIndex = new(arrayMax(vt->sites), int) ;
  DICT *gtDict = dictCreate (4*arrayMax(vt->sites)) ;
  HASH *posHash = hashCreate (4*arrayMax(vt->sites)) ;
  BOOL *isAnc1 = new (arrayMax(vt->sites), BOOL) ;
  double *llSite = new (arrayMax(vt->sites), double) ;

  Array hMiss = arrayCreate (256, int) ; // histogram
  int nWithMiss = 0, totMiss = 0 ;
  for (i = 0 ; i < arrayMax(vt->sites) ; ++i)
    { VcfSite *s = arrp(vt->sites, i, VcfSite) ;
      if (!hashAdd (posHash, HASH_INT(POS_HASH(s)), &k) || k != i)
	die ("problem hashing position %d", i) ;
	  
      int nMiss = 0, n0 = 0, n1 = 0, nBad = 0 ;
      for (j = 0 ; j < arrayMax (vt->samples) ; ++j)
	switch (s->gt[j])
	  {
	  case 1: ++nMiss ; break ; case 2: ++n0 ; break ; case 3: ++n1 ; break ;
	  default: ++nBad ;
	  }
      if (nBad) { fprintf (stderr, "bad gt %d for site %d\n", nBad, i) ; continue ; }
      if (!n1) { fprintf (stderr, "no gt1 for site %d\n", i) ; continue ; }
      if (nMiss) { ++nWithMiss ; totMiss += nMiss ; }

      if (dictAdd (gtDict, s->gt, &siteIndex[i])) // a new site pattern
	{ char R = toupper(s->ref), A = toupper(s->alt) ;
	  if ((R == 'A' && A == 'G') || (R == 'C' && A == 'T') ||
	      (R == 'G' && A == 'A') || (R == 'T' && A == 'C'))
	    scores[siteIndex[i]] = treeBuildScores (t, s->gt, tree2vcf, transitionEdges) ;
	  else
	    scores[siteIndex[i]] = treeBuildScores (t, s->gt, tree2vcf, transversionEdges) ;
	  if (isVerbose) ++array(hMiss, nMiss, int) ;
	}
      TreeScore *score = scores[siteIndex[i]] ;
      if (score[0].below.s0 < score[0].below.s1)
	{ isAnc1[i] = TRUE ;
	  llSite[i] = score[0].below.s1 + log(1. + exp(score[0].below.s0-score[0].below.s1)) ;
	}
      else
	{ isAnc1[i] = FALSE ;
	  llSite[i] = score[0].below.s0 + log(1. + exp(score[0].below.s1-score[0].below.s0)) ;
	}
    }
  printf ("built scores for %d gt patterns for %d sites\n", dictMax(gtDict), arrayMax(vt->sites)) ;
  if (nWithMiss) printf ("%d sites had missing genotypes, mean %.1f\n",
			 nWithMiss, totMiss / (double) nWithMiss) ;
  if (isVerbose)
    for (i = 0 ; i < arrayMax (hMiss) ; ++i)
      if ((j = arr(hMiss,i,int))) printf ("  %d site patterns with %d genotypes missing\n", j, i) ;
  arrayDestroy (hMiss) ;
  timeUpdate (stdout) ;

  if (fB)
    { for (i = 0 ; i < arrayMax(vt->sites) ; ++i)
	{ VcfSite *s = arrp(vt->sites, i, VcfSite) ;
	  double best = 0., best2 = 0. ; int kBest, kBest2 ;
	  TreeScore *score = scores[siteIndex[i]] ;
	  for (k = 1 ; k < arrayMax(t->a) ; ++k)
	    { double x = isAnc1[i] ? score[k].below.s0 + score[k].above.s1
		                   : score[k].below.s1 + score[k].above.s0 ;
	      if (!best || x > best) { best2 = best ; best = x ; kBest = k ; }
	      else if (!best2 || x > best2) { best2 = x ; kBest2 = k ; }
	    }
	  fprintf (fB, "%s\t%d\t%c\t%c\t%d\t%d\t%8.2f\t%d\t%8.2f\n",
		   vt->seqName, s->pos, s->ref, s->alt, isAnc1[i], kBest, best, kBest2, best2-best) ;
	}
      fclose (fB) ;
      printf ("written branch assignments for %d sites to file\n", arrayMax(vt->sites)) ;
      timeUpdate (stdout) ;
    }

  if (vq)
    { if (strcmp (vq->seqName, vt->seqName))
	die ("query seqName %s != reference %s", vq->seqName, vt->seqName) ;

      printf ("query %d samples at %d sites to tree\n", dictMax(vq->samples), arrayMax(vq->sites)) ;
      
      int *siteMap = new (arrayMax(vq->sites), int) ; // index of vq site in vt site list
      int nSitesMapped = arrayMax(vq->sites) ;
      for (i = 0 ; i < arrayMax (vq->sites) ; ++i)
	{ VcfSite *s = arrp (vq->sites, i, VcfSite) ;
	  if (hashFind (posHash, HASH_INT(POS_HASH(s)), &k))
	    siteMap[i] = k ;
	  else
	    { siteMap[i] = -1 ; --nSitesMapped ; }
	}
      printf ("  matched %d sites\n", nSitesMapped) ;
      if (!nSitesMapped) die ("no shared sites to map query samples with\n") ;

      double *qScore = new (arrayMax(t->a), double) ;
      for (j = 0 ; j < dictMax (vq->samples) ; ++j)
	{ bzero (qScore, arrayMax(t->a)*sizeof(double)) ;
	  double baseScore = 0. ;
	  int n = 0 ;
	  for (i = 0 ; i < arrayMax (vq->sites) ; ++i)
	    if (siteMap[i] >= 0)
	      { TreeScore *score = scores[siteIndex[siteMap[i]]] ;
		VcfSite *sq = arrp (vq->sites, i, VcfSite) ;
		if (sq->gt[j] > 1)
		  { ++n ;
		    if (sq->gt[j] == 2) // ref
		      for (k = 0 ; k < arrayMax(t->a) ; ++k)
			qScore[k] += score[k].below.s0 + score[k].above.s0 ;
		    else if (sq->gt[j] == 3) // alt
		      for (k = 0 ; k < arrayMax(t->a) ; ++k)
			qScore[k] += score[k].below.s1 + score[k].above.s1 ;
		    baseScore += llSite[siteMap[i]] ;
		  }
	      }
	  if (!n)
	    { printf ("  %s no data to map\n", dictName(vq->samples,j)) ;
	      continue ;
	    }
	  double best = 0. ; int kBest ;
	  for (k = 0 ; k < arrayMax(t->a) ; ++k)
	    if (!best || qScore[k] > best) { best = qScore[k] ; kBest = k ; }
	  printf ("  %s\tn %d\tbranch %d\tscore %.2f %.2f\n", dictName(vq->samples,j), n,
		  kBest, best - baseScore, (best - baseScore)/n) ;
	  if (isVerbose)
	    { double diff = log(100.0) ;
	      for (k = 0 ; k < arrayMax(t->a) ; ++k)
		if (k != kBest && qScore[k] > best-diff)
		  printf ("    sub-optimal branch %d\trelative-score %.3f\n",
			  k, exp(qScore[k] - best)) ;
	    }
	}
      vcfDestroy (vq) ;
      free (qScore) ;
      timeUpdate (stdout) ;
    }

  treeDestroy (t) ;
  vcfDestroy (vt) ;
  dictDestroy (gtDict) ;
  hashDestroy (posHash) ;

  timeTotal (stdout) ;
  
  return 0 ;
}

/************ end ************/
