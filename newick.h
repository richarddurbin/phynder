/*  File: newick.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 09:59 2019 (rd109)
 * Created: Fri Nov 15 23:58:04 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "tree.h"

TreeNode *readBinaryNewickTree (FILE *f) ;
void writeBinaryNewickTree (FILE *f, TreeNode *n) ;

/******************************************************/
