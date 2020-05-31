/*  File: treevcf.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: based on pbwtHtslib.c
 * Exported functions:
 * HISTORY:
 * Last edited: May 14 03:45 2020 (rd109)
 * Created: Sun Nov 17 19:41:19 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vcf.h"
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
#include <ctype.h>		/* for toupper() */

Vcf *vcfRead (char *filename, int *multi, int *nonSNP)  /* read GTs from vcf/bcf using htslib */
{
  int i, j, k ;

  bcf_srs_t *sr = bcf_sr_init () ;
  if (!bcf_sr_add_reader (sr, filename)) die ("failed to open vcf file", filename) ;

  Vcf *v = new0 (1, Vcf) ;
  v->samples = dictCreate (1024) ;
  v->sites = arrayCreate (1024, VcfSite) ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  int nSample = bcf_hdr_nsamples(hr) ;
  for (i = 0 ; i < nSample ; ++i)
    { if (!dictAdd (v->samples, hr->samples[i], &k))
	fprintf (stderr, "repeat sample name %s in vcf file\n", hr->samples[i]) ;
      assert (k == i) ;
    }
  assert (nSample == dictMax(v->samples)) ;

  int ngt_arr = 0, *gt_arr = NULL ;
  *multi = 0, *nonSNP = 0 ;
  while (bcf_sr_next_line (sr)) 
    { bcf1_t *line = bcf_sr_get_line(sr,0) ;
      if (!v->seqName)
	v->seqName = strdup (bcf_seqname(hr,line)) ;
      else if (strcmp (v->seqName, bcf_seqname(hr,line)))
	die ("multiple sequence names in VCF file %s != %s", v->seqName, bcf_seqname(hr,line)) ;
      if (line->n_allele == 1) continue ; // ignore lines with no alts
      if (line->n_allele > 2) { ++*multi ; continue ; }
      if (line->d.allele[0][1] || line->d.allele[1][1]) { ++*nonSNP ; continue ; }
      
      // get GTs
      int ngt = bcf_get_genotypes(hr, line, &gt_arr, &ngt_arr) ;
      if (ngt <= 0) continue ;  // -1 is used if GT is not in the FORMAT
      if (ngt != nSample)
	die ("wrong number of genotypes %d at VCF %s %lld", ngt, v->seqName, line->pos+1) ;

      VcfSite *s = arrayp(v->sites,arrayMax(v->sites),VcfSite) ;
      s->pos = line->pos + 1 ;       // bcf coordinates are 0-based
      s->ref = tolower(*line->d.allele[0]) ;
      s->alt = tolower(*line->d.allele[1]) ;
      s->gt = new (nSample+1, char) ; // zero-terminate so can use as a string in dict package
      for (i = 0 ; i < nSample ; i++)
	{ if (gt_arr[i] == bcf_int32_vector_end) 
	    die ("unexpected end of genotype vector in VCF") ;
	  if (gt_arr[i] == bcf_gt_missing)
	    s->gt[i] = 1 ;
	  else
	    s->gt[i] = 2 + bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 2 or 3
	}
      s->gt[nSample] = 0 ;
    }

  bcf_sr_destroy (sr) ;
  
  return v ;
}

void vcfDestroy (Vcf *v)
{
  dictDestroy (v->samples) ;
  int i ;
  for (i = 0 ; i < arrayMax(v->sites) ; ++i) free (arrp(v->sites,i,VcfSite)->gt) ;
  arrayDestroy (v->sites) ;
  free (v->seqName) ;
  free (v) ;
}

/******* end of file ********/
