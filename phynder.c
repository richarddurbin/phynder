/*  File: phynder.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: likelihood applications on trees for ancient Y data and more
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 13 11:44 2020 (rd109)
 * Created: Sun Nov 17 19:52:20 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "tree.h"
#include "newick.h"
#include "vcf.h"
#include "ONElib.h"
#include <ctype.h>
#include <math.h>
#include <stdarg.h>

char *schemaText =
  "1 3 def 1 0  schema for phynder\n"
  ".\n"
  "P 3 var                    variant file\n"
  "S 3 snp                    SNP file\n"
  "D c 2 3 INT 6 STRING          chromosome\n"
  "D V 2 3 INT 11 STRING_LIST    pos, ref string, alt string\n"
  "D A 1 4 CHAR                  0 for ref ancestral, 1 for alt ancestral\n"
  "D L 1 4 REAL                  log likelihood of the site\n"
  "D B 2 3 INT 4 REAL            best branch assignment, score\n"
  "D S 2 3 INT 4 REAL            secondary branch assignment, score\n"
  ".\n"
  "P 3 smp                    sample file\n"
  "D I 1 6 STRING                sample name\n"
  "D N 1 3 INT                   number of sites\n"
  "D B 3 3 INT 4 REAL 4 REAL     best branch, posterior, score\n"
  "D S 3 3 INT 4 REAL 4 REAL     suboptimal branch, posterior, score\n"
  "D C 1 3 INT                   clade\n"
  ".\n"
  "P 3 nul                    empty file - comments only\n" ;

static double transitionRate = 1.33 ; // expect ratio 4:1 for ts:tv rates
static double transversionRate = 0.33 ; // these values give mean 1.0 for 2:1 ts:tv events
static double siteThreshold = -10.0 ; // log likelihood threshold to accept a site pattern
static BOOL   isUltrametric = FALSE ;
static BOOL   isVerbose = FALSE ;
static int    calcMode = 0 ; // calculation mode
static double posteriorThreshold = 0. ;

static int baseMap[128] ;
#define POS_HASH(s) (((s)->pos << 4) + (baseMap[(s)->ref] << 2) + baseMap[(s)->alt])

static void usage (void)
{
  fprintf (stderr, "usage: phynder [options] <newick tree> <vcf for tree>\n") ;
  fprintf (stderr, "  -o <output filename>     ONEfile name; default is stdout [-]\n") ;
  fprintf (stderr, "  -b                       write binary (can read with ONEview)\n") ;
  fprintf (stderr, "  -q <query vcf>           output query assignments\n") ;
  fprintf (stderr, "  -p <posterior threshold> print out suboptimal branches and clade [%.0f]\n", posteriorThreshold) ;
  fprintf (stderr, "  -B                       output branch positions of tree variants\n") ;
  fprintf (stderr, "  -T <thresh>              site likelihood threshold, zero means no threshold [%.1f]\n", siteThreshold) ;
  fprintf (stderr, "  -ts <transition rate>    [%.4f]\n", transitionRate) ;
  fprintf (stderr, "  -tv <transversion rate>  [%.4f]\n", transversionRate) ;
  fprintf (stderr, "  -C <calc_mode>           [%d] calculation mode\n", calcMode) ;
  fprintf (stderr, "                           calc_mode 0: LL both ends of edge match\n") ;
  fprintf (stderr, "                           calc_mode 1: -LL both ends of edge mismatch\n") ;
  fprintf (stderr, "  -U                       make tree ultrametric - all leaves equidistant from root\n") ;
  fprintf (stderr, "  -v                       verbose - print extra info\n") ;
  fprintf (stderr, "  -h                       print this message\n") ;
  exit (0) ;
}

static char *commandLine (int argc, char **argv)
{
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  char *buf = new (totLen, char) ;
  strcpy (buf, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (buf, " ") ; strcat (buf, argv[i]) ; }
  return buf ;
}

void onePrintf (OneFile *vf, char *format, ...)
{
  static char buf[1024] ;
  va_list args ;

  oneWriteLine (vf, '.', 0, 0) ;

  va_start (args, format) ;
  vsprintf (buf, format, args) ;
  va_end (args) ;

  oneWriteComment (vf, buf) ;
}

int main (int argc, char *argv[])
{
  char      *queryFile = 0, *outFile = "-" ;
  bool       isBranchOut = false, isBinary = false ;
 
  timeUpdate (0) ;

  baseMap['a'] = baseMap['A'] = 0 ;
  baseMap['c'] = baseMap['C'] = 1 ;
  baseMap['g'] = baseMap['G'] = 2 ;
  baseMap['t'] = baseMap['T'] = 3 ;
  baseMap['0'] = 0 ;
  baseMap['1'] = 1 ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ; // absorb executable name
  while (argc)
    if (!strcmp (*argv, "-o") && argc > 1)
      { outFile = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-b"))
      { isBinary = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-q") && argc > 1)
      { queryFile = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-B"))
      { isBranchOut = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-ts") && argc > 1)
      { transitionRate = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-tv") && argc > 1)
      { transversionRate = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-T") && argc > 1)
      { siteThreshold = atof (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-p") && argc > 1)
      { posteriorThreshold = atof (argv[1]) ; argc -= 2 ; argv += 2 ;
	if (posteriorThreshold < 0. || posteriorThreshold >= 1.)
	  die ("posterior %f must be between 0 and 1", posteriorThreshold) ;
      }
    else if (!strcmp (*argv, "-v"))
      { isVerbose = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-C") && argc > 1)
      { calcMode = atoi (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-h"))
      usage () ;
    else if (argc != 2)
      die ("command line error at %s - run without arguments for usage", *argv) ;
    else
      break ;

  if (!argc) usage () ;
  
  OneSchema *vs = oneSchemaCreateFromText (schemaText) ;
  OneFile *vf ; 
  if (isBranchOut) vf = oneFileOpenWriteNew (outFile, vs, "snp", isBinary, 1) ;
  else if (queryFile) vf = oneFileOpenWriteNew (outFile, vs, "smp", isBinary, 1) ;
  else  vf = oneFileOpenWriteNew (outFile, vs, "nul", isBinary, 1) ;
  oneSchemaDestroy (vs) ;
  if (!vf) die ("failed to open output file %s", outFile) ;
  oneAddProvenance (vf, "phynder", "1.0", command, 0) ;
  oneWriteHeader (vf) ;

  FILE *f ;
  if (!(f = fopen (argv[0], "r"))) die ("failed to open newick tree file %s", argv[0]) ;
  TreeNode *n = readBinaryNewickTree (f) ;
  fclose (f) ;
  Tree *t = treeCreate (n) ;
  treeNodeDestroy (n) ;
  onePrintf (vf, "read tree with %d nodes and %d leaves",
	       arrayMax(t->a), (arrayMax(t->a)+1)/2) ;
  fprintf (stderr, "read tree with %d nodes: ", arrayMax(t->a)) ;
  timeUpdate (stderr) ;

  int multi = 0, nonSNP = 0 ;
  Vcf *vt = vcfRead (argv[1], &multi, &nonSNP) ;
  onePrintf (vf, "read %d sites for %d samples in tree",
	     arrayMax(vt->sites), dictMax(vt->samples)) ;
  if (multi || nonSNP)
    onePrintf (vf, "  ignored %d multiple-allelele and %d non-SNP sites",
	       multi, nonSNP) ;
  fprintf (stderr, "read vcf for tree at %d sites: ", arrayMax(vt->sites)) ;
  timeUpdate (stderr) ;

  if (isUltrametric)
    { treeBalance (t) ;
      fprintf (stderr, "balanced tree: ") ;
      timeUpdate (stderr) ;
    }

  double worstTransition, worstTransversion ;
  Array transitionEdges = treeBuildEdges (t, transitionRate, &worstTransition) ;
  Array transversionEdges = treeBuildEdges (t, transversionRate, &worstTransversion) ;
  if (isVerbose)
    onePrintf (vf, "  worst transition transversion %8.2f %8.2f",
		 worstTransition, worstTransversion) ;

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
  BOOL isMissingSample = FALSE ;
  for (i = 0 ; i < arrayMax(t->a) ; ++i)
    if (!arrp(t->a,i,TreeElement)->left && !leafFound[i])
      { char *name = dictName (t->nameDict, i) ;
	if (!dictFind (vt->samples, name, &j))
	  { fprintf (stderr, "missing sample '%s' in tree vcf\n", name) ;
	    isMissingSample = TRUE ;
	  }
	else die ("problem connecting tree leaf %s i %d to vcf sample j %d", name, i, j) ;
      }
  if (isMissingSample) die ("missing samples in tree vcf file") ;
  free (leafFound) ;
  
  // build scores for each site for each node in tree
  // indirection via siteIndex so that info for sites with identical genotypes is shared
  // NB if -B to output branches then can't also fit queries
  if (isBranchOut) calcMode = -1 ;
  LogLikelihood **scores = new0(arrayMax(vt->sites), LogLikelihood*) ; // more than needed but cheap
  int *siteIndex = new0(arrayMax(vt->sites), int) ;
  DICT *gtDict = dictCreate (4*arrayMax(vt->sites)) ;
  dictAdd (gtDict, "burn the first siteIndex entry", 0) ; // needed so siteIndex[i] = 0 for bad sites
  HASH *posHash = hashCreate (4*arrayMax(vt->sites)) ;
  BOOL *isAnc1 = new (arrayMax(vt->sites), BOOL) ;
  double *llSite = new (arrayMax(vt->sites), double) ;
  BOOL *isTransition = new (arrayMax(vt->sites), BOOL) ;

  Array hMiss = arrayCreate (256, int) ; // histogram
  int nWithMiss = 0, totMiss = 0, nThresh = 0, nBadGT = 0, nMonomorphic = 0, nGood = 0 ;
  for (i = 0 ; i < arrayMax(vt->sites) ; ++i)
    { VcfSite *s = arrp(vt->sites, i, VcfSite) ;
      if (!hashAdd (posHash, HASH_INT(POS_HASH(s)), &k) || k != i)
	die ("duplicated site in VCF: %d %c %c", s->pos, s->ref, s->alt) ;

      { char R = toupper(s->ref), A = toupper(s->alt) ;
	isTransition[i] = ((R == 'A' && A == 'G') || (R == 'C' && A == 'T') ||
			   (R == 'G' && A == 'A') || (R == 'T' && A == 'C')) ;
      }	
      
      int nMiss = 0, n0 = 0, n1 = 0, nBad = 0 ;
      for (j = 0 ; j < arrayMax (vt->samples) ; ++j)
	switch (s->gt[j])
	  {
	  case 1: ++nMiss ; break ; case 2: ++n0 ; break ; case 3: ++n1 ; break ;
	  default: ++nBad ;
	  }
      if (nBad)
	{ if (isVerbose)
	    fprintf (stderr, "%d bad genotypes for site %d - drop site\n", nBad, i) ;
	  ++nBadGT ;
	  continue ;
	}
      if (!n0 || !n1)
	{ if (isVerbose)
	    fprintf (stderr, "site %d is monomorphic gt0 %d gt1 %d\n", i, n0, n1) ;
	  ++nMonomorphic ;
	  continue ;
	}
      if (nMiss) { ++nWithMiss ; totMiss += nMiss ; }

      if (dictAdd (gtDict, s->gt, &k)) // a new site pattern
	{ if (isTransition[i])
	    scores[k] = treeBuildScores (t, s->gt, tree2vcf, transitionEdges, calcMode) ;
	  else
	    scores[k] = treeBuildScores (t, s->gt, tree2vcf, transversionEdges, calcMode) ;
	  if (isVerbose) ++array(hMiss, nMiss, int) ;
	}
      LogLikelihood *score = scores[k] ;
      if (score[0].s0 < score[0].s1)
	{ isAnc1[i] = TRUE ;
	  llSite[k] = score[0].s1 + log(1. + exp(score[0].s0-score[0].s1)) ;
	}
      else
	{ isAnc1[i] = FALSE ;
	  llSite[k] = score[0].s0 + log(1. + exp(score[0].s1-score[0].s0)) ;
	}
      double thisThresh ;
      if (isTransition[i]) thisThresh = worstTransition ; else thisThresh = worstTransversion ;
      if (siteThreshold && llSite[k] < siteThreshold + thisThresh)
	{ ++nThresh ;
	  continue ;
	}

      siteIndex[i] = k ;  // only set the siteIndex for good sites, else 0 indicating bad
      ++nGood ;
    }
  onePrintf (vf, "built scores for %d gt patterns", dictMax(gtDict)) ;
  onePrintf (vf, "using %d good sites", nGood) ;
  if (nMonomorphic) onePrintf (vf, "  %d sites rejected because monomorphic", nMonomorphic) ;
  if (nBadGT) onePrintf (vf, "  %d sites rejected because they had bad genotpes", nBadGT) ;
  if (nThresh) onePrintf (vf, "  %d sites rejected because likelihood below threshold %.1f",
			  nThresh, siteThreshold) ;
  if (nWithMiss) onePrintf (vf, "  %d sites had missing genotypes, mean %.1f",
			    nWithMiss, totMiss / (double) nWithMiss) ;
  if (isVerbose)
    for (i = 0 ; i < arrayMax (hMiss) ; ++i)
      if ((j = arr(hMiss,i,int)))
	fprintf (stderr, "  %d site patterns with %d genotypes missing\n", j, i) ;
  arrayDestroy (hMiss) ;
  fprintf (stderr, "built score table: ") ;
  timeUpdate (stderr) ;

  double *logPrior = new0 (arrayMax(t->a), double) ;

  // code for finding the most likely branch(es) for mutations per site
  if (isBranchOut)
    { onePrintf (vf, "") ; // end of header reports
      oneInt(vf,0) = arrayMax(vt->sites) ;
      oneWriteLine (vf, 'c', strlen(vt->seqName), vt->seqName) ;

      char *buf = new0 (4, char) ;
      for (i = 0 ; i < arrayMax(vt->sites) ; ++i)
	{ if (!siteIndex[i]) continue ;
	  VcfSite *s = arrp(vt->sites, i, VcfSite) ;
	  double best = 0., best2 = 0. ; int kBest, kBest2 ;
	  LogLikelihood *score = scores[siteIndex[i]] ;
	  for (k = 1 ; k < arrayMax(t->a) ; ++k)
	    { double x = isAnc1[i] ? score[k].s1 : score[k].s0 ;
	      if (!best || x > best) { best2 = best ; kBest2 = kBest ; best = x ; kBest = k ; }
	      else if (!best2 || x > best2) { best2 = x ; kBest2 = k ; }
	    }
	  oneInt(vf,0)=s->pos ; buf[0]=s->ref ; buf[2]=s->alt ; oneWriteLine (vf, 'V', 2, buf) ;
	  oneChar(vf,0) = isAnc1[i] ? '1' : '0' ; oneWriteLine (vf, 'A', 0, 0) ;
	  oneReal(vf,0) = llSite[siteIndex[i]] ; oneWriteLine (vf, 'L', 0, 0) ;
	  oneInt(vf,0) = kBest ; oneReal(vf,1) = best ; oneWriteLine (vf, 'B', 0, 0) ;
	  oneInt(vf,0) = kBest2 ; oneReal(vf,1) = best2 ; oneWriteLine (vf, 'S', 0, 0) ;
	}

      fprintf (stderr, "assigned %d branches: ", arrayMax(vt->sites)) ;
      timeUpdate (stderr) ;
    }
  
  // code for finding the most likely branch(es) for new samples
  else if (queryFile)
    { Vcf *vq = vcfRead (queryFile, &multi, &nonSNP) ;
      onePrintf (vf, "read %d sites for %d samples in query",
		 arrayMax(vq->sites), dictMax(vq->samples)) ;
      if (multi || nonSNP)
	onePrintf (vf, "  ignored %d multiple-allelele and %d non-SNP sites",
		   multi, nonSNP) ;
      fflush (vf->f) ;
      onePrintf (vf, "") ; // end of header reports
  
      if (strcmp (vq->seqName, vt->seqName))
	die ("query seqName %s != reference %s", vq->seqName, vt->seqName) ;

      int *siteMap = new0 (arrayMax(vq->sites), int) ; // index of vq site in vt site list
      int nSitesMapped = 0 ;
      for (i = 0 ; i < arrayMax (vq->sites) ; ++i)
	{ VcfSite *s = arrp (vq->sites, i, VcfSite) ;
	  if (hashFind (posHash, HASH_INT(POS_HASH(s)), &k) && siteIndex[k])
	    { siteMap[i] = siteIndex[k] ;
	      ++nSitesMapped ;
	    }
	  else if (isVerbose)
	    fprintf (stderr, "  can't find site in tree: %d %c %c\n",s->pos, s->ref, s->alt) ;
	}
      if (!nSitesMapped) die ("no shared sites with which to map query samples") ;

  //      printf ("file contains %d samples with data at %d sites,", 
  //	      dictMax(vq->samples), arrayMax(vq->sites)) ;
  //      printf (" matching %d sites in the tree\n", nSitesMapped) ;

      double *qScore = new (arrayMax(t->a), double) ;
      double *qCladeTotal = new (arrayMax(t->a), double) ;
      int nMapped = 0 ;
      for (j = 0 ; j < dictMax (vq->samples) ; ++j)
	{ memcpy (qScore, logPrior, arrayMax(t->a)*sizeof(double)) ;
	  double baseScore = 0. ;
	  int nSite = 0 ;
	  for (i = 0 ; i < arrayMax (vq->sites) ; ++i)
	    if (siteMap[i])
	      { LogLikelihood *score = scores[siteMap[i]] ;
		VcfSite *sq = arrp (vq->sites, i, VcfSite) ;
		if (sq->gt[j] > 1)
		  { ++nSite ;
		    // around 191231 started writing the line below
		    // treeQueryAddSite (qScore, sq->gt[j]-2, scores[siteMap[i]],
		    if (sq->gt[j] == 2) // ref
		      for (k = 1 ; k < arrayMax(t->a) ; ++k)
			qScore[k] += score[k].s0 ;
		    else if (sq->gt[j] == 3) // alt
		      for (k = 1 ; k < arrayMax(t->a) ; ++k)
			qScore[k] += score[k].s1 ;
		    baseScore += llSite[siteMap[i]] ;
		  }
	      }
	  if (!nSite)
	    { fprintf (stderr, "query %s no data to map\n", dictName(vq->samples,j)) ;
	      continue ;
	    }
	  ++nMapped ;
	  double best = 0. ; int kBest ;
	  for (k = 1 ; k < arrayMax(t->a) ; ++k)
	    if (!best || qScore[k] > best) { best = qScore[k] ; kBest = k ; }
	  double total = 0. ;
	  for (k = 1 ; k < arrayMax(t->a) ; ++k)
	    if (qScore[k] > best-15.) total += exp(qScore[k] - best) ;

	  oneWriteLine (vf, 'I', strlen(dictName(vq->samples, j)), dictName(vq->samples, j)) ;
	  oneInt(vf,0) = nSite ; oneWriteLine (vf, 'N', 0, 0) ;
	  oneInt(vf,0) = kBest ; oneReal(vf,1) = 1./total ;
	    oneReal(vf,2) = best - baseScore ; oneWriteLine (vf, 'B', 0, 0) ;
	  
	  if (posteriorThreshold)
	    { double diff = log(posteriorThreshold) ;
	      for (k = 1 ; k < arrayMax(t->a) ; ++k)
		if (k != kBest && qScore[k] > best+diff)
		  { oneInt(vf,0) = k ; oneReal(vf,1) = exp(qScore[k] - best) / total ;
		      oneReal(vf,2) = qScore[k] - baseScore ; oneWriteLine (vf, 'S', 0, 0) ;
		  }
	      double cladeThresh = total * (1-posteriorThreshold) ;
	      for (k = arrayMax(t->a) ; k-- ; ) // need reverse order for post-order
		{ TreeElement *e = arrp(t->a, k, TreeElement) ;
		  if (e->left)
		    qCladeTotal[k] = qCladeTotal[e->left] + qCladeTotal[e->right] ;
		  else
		    qCladeTotal[k] = 0. ;
		  if (qScore[k] > best-15.)
		    qCladeTotal[k] += exp(qScore[k] - best) ;
		  if (qCladeTotal[k] > cladeThresh) break ;
		}
	      oneInt(vf,0) = k ; oneWriteLine (vf, 'C', 0, 0) ;
	    }
	}
      vcfDestroy (vq) ;
      free (qScore) ;
      free (qCladeTotal) ;

      fprintf (stderr, "mapped %d samples: ", nMapped) ;
      timeUpdate (stderr) ;
    }

  oneFileClose (vf) ;
  
  treeDestroy (t) ;
  vcfDestroy (vt) ;
  dictDestroy (gtDict) ;
  hashDestroy (posHash) ;

  fprintf (stderr, "total resource usage: ") ;
  timeTotal (stderr) ;
  
  return 0 ;
}

/************ end ************/
