/*  File: vcf.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 19 18:28 2019 (rd109)
 * Created: Sun Nov 17 23:25:45 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef struct {
  char *gt ;			// 1 = miss, 2 = 0, 3 = 1, 0-terminated: length dictMax(samples)+1
  int pos ;
  char ref, alt ;		// lower case - only parse binary SNPs for now
} VcfSite ;

typedef struct {
  DICT *samples ;
  Array sites ;			// of VcfSite
  char *seqName ;
} Vcf ;

Vcf *vcfRead (char *filename) ;
void vcfDestroy (Vcf *v) ;

/****************************************/
