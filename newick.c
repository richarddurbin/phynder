/*  File: newick.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: based on Matt Rasmussen's http://compbio.mit.edu/spimap/pub/spimap/src/newick.cpp
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 25 10:48 2020 (rd109)
 * Created: Fri Nov 15 19:17:14 2019 (rd109)
 *-------------------------------------------------------------------
 */

/* from https://en.wikipedia.org/wiki/Newick_format

   The grammar nodes
     Tree: The full input Newick Format for a single tree
     Subtree: an internal node (and its descendants) or a leaf node
     Leaf: a node with no descendants
     Internal: a node and its one or more descendants
     BranchSet: a set of one or more Branches
     Branch: a tree edge and its descendant subtree.
     Name: the name of a node
     Length: the length of a tree edge.

   The grammar rules
     Tree → Subtree ";" | Branch ";"
     Subtree → Leaf | Internal
     Leaf → Name
     Internal → "(" BranchSet ")" Name
     BranchSet → Branch | Branch "," BranchSet
     Branch → Subtree Length
     Name → empty | string
     Length → empty | ":" number

   Whitespace (spaces, tabs, carriage returns, and linefeeds) within number is prohibited. 
   Whitespace within string is often prohibited.
   Whitespace elsewhere is ignored. 
   The punctuation characters from the grammar (semicolon, parentheses, comma, and colon) are prohibited in string.
   The Tree --> Branch ";" production makes the entire tree descendant from nowhere, which can be nonsensical, and is sometimes prohibited.

   I allow whitespace in string, but 
*/

#include "newick.h"
#include <ctype.h>

static BOOL matchChar (FILE *f, char t)
{
  char c ;
  do { c = getc(f) ; } while (isspace (c)) ;
  if (feof(f)) die ("unexpected end of newick file") ;
  if (c == t) return TRUE ;
  ungetc (c, f) ; return FALSE ;
}

static TreeNode *readBinarySubTree (FILE *f, TreeNode *parent)
{
  static char buf[1024] ; // NB limit to node name size
  TreeNode *n = new0 (1, TreeNode) ;
  n->parent = parent ;
  if (matchChar (f, '('))
    { n->left = readBinarySubTree (f, n) ;
      if (!matchChar (f, ',')) die ("missing comma in newick file") ;
      n->right = readBinarySubTree (f, n) ;
      if (!parent && matchChar(f, ',')) n->parent = readBinarySubTree (f, n) ;
      if (!matchChar (f, ')')) die ("missing close parenthesis in newick file") ;
    }
  while (TRUE)
    { matchChar (f, 0) ;  // absorb whitespace and guarantee next char
      char c = fgetc (f) ;
      switch (c)
	{
	case ',': case ')': case ';': // end of node
	  ungetc (c, f) ; // return separator to file
	  return n ;
	case ':': // parse length
	  if (n->length) die ("double length specification in newick file") ;
	  fscanf (f, "%lf", &n->length) ;
	  if (!n->length)
	    { fprintf (stderr, "branch length 0 converted to 1e-8\n") ;
	      n->length = 1.0e-8 ;
	    }
	  break ;
	default:
	  if (n->name || n->length) die ("incorrect name syntax in newick file") ;
	  ungetc (c, f) ; // return the char which starts the name
	  fscanf (f, "%[^,;:()]", buf) ;
	  { char *c = buf ; while (*c) ++c ;
	    --c ; while (isspace (*c)) --c ; *++c = 0 ; // trim trailing whitespace
	    n->name = new (c-buf+1, char) ; strcpy (n->name, buf) ;
	  }
	}
    }
  return 0 ; // actually the code can't reach here, but have to keep dumb compilers happy
}

TreeNode *readBinaryNewickTree (FILE *f)
{
  TreeNode *n = readBinarySubTree (f, 0) ;
  n->isRoot = TRUE ;
  if (!matchChar (f, ';')) die ("newick tree parse failed to end in ';'") ;
  return n ;
}

/************************/

static void writeBinarySubTree (FILE *f, TreeNode *n, int depth)
{
  int i ;
  if (n->left)
    { ++depth ;
      fputc ('(', f) ; writeBinarySubTree (f, n->left, depth) ;
      fputc (',', f) ; fputc ('\n', f) ; 
      for (i = depth ; i-- ;) fputc (' ', f) ;
      writeBinarySubTree (f, n->right, depth) ;
      if (!depth && n->parent)
	{ fputc (',', f) ; fputc ('\n', f) ; 
	  for (i = depth ; i-- ;) fputc (' ', f) ;
	  writeBinarySubTree (f, n->parent, depth) ;
	}
      --depth ;
      fputc ('\n', f) ;
      for (i = depth ; i-- ;) fputc (' ', f) ;
      fputc (')', f) ;
    }
  if (n->name) fprintf (f, " %s", n->name) ;
  if (n->length) fprintf (f, " :%.12g", n->length) ;
}

void writeBinaryNewickTree (FILE *f, TreeNode *n)
{ writeBinarySubTree (f, n, 0) ;
  fprintf (f, ";\n") ;
}

/************************/

#ifdef TEST

char *testString = "(HG02982:0.016947244492149634,((HG02613:0.0010573526625989628,(HG02666:4.723421115815363E-5,HG02645:1.742760356723723E-4):0.0011135732100915653):0.012120421043158552,(((ERS474081:0.001336208244420229,LP6005441-DNA_B11:0.001370666848282548):0.01170993621272617,(SS6004473:0.005789981089206839,(SS6004480:4.4150145785174996E-4,(LP6005443-DNA_B09:2.6911902024052275E-4,ERS474093:2.5941043238843567E-4):8.638445992146734E-5):0.005507103801356057):0.006720527891845726):3.161207493461403E-4,((ERS474111:0.008656365429731316,(((ERS474074:2.4380154420359373E-5,ERS474068:1.2261853507687659E-5):0.00398677000185161,(GS000035246-ASM:0.002220414190312936,(ERS474839:0.0010872828548510946,(ERS474089:3.904595888214006E-4,LP6005441-DNA_G02:3.6611754204254887E-4):5.897901564099368E-4):0.001173473446796084):0.0015668311391512264):0.0019645650955184175,((ERS474059:1.362675602461401E-4,ERS474043:4.867408542482202E-5):0.0054770977743243965,((LP6005441-DNA_A08:0.0038238205621545333,(LP6005443-DNA_G08:8.564390907284461E-5,LP6005441-DNA_A11:7.315583889601351E-5):0.0042030713667862325):3.969609422279953E-4,(ERS474080:0.004149933502997694,((ERS474070:1.2274972773330256E-5,ERS474076:1.2280367044419094E-5):0.002982026277018126,(LP6005441-DNA_H02:0.002345644645735539,(ERS474042:1.589530822464838E-4,(ERS474063:1.2207724657171265E-5,ERS474058:6.101470470837658E-5):4.864989842363193E-5):0.0024316447192818703):5.120975621530424E-4):0.001131454945690235):2.3055943918854488E-4):7.985375925601626E-4):9.600279166768685E-4):0.0029768094443579085):2.793373660597148E-4,(((('/lustre/scratch119/humgen/teams/tyler-smith/users/ph8/BigTree/D0/remap/yhg-de7119651_chrY_hg38full-woD_srt.bam':1.2286034332632118E-4,('/lustre/scratch119/humgen/teams/tyler-smith/users/ph8/BigTree/D0/remap/yhg-de7119652_chrY_hg38full-woD_srt.bam':1.2302160666773088E-4,'/lustre/scratch119/humgen/teams/tyler-smith/users/ph8/BigTree/D0/remap/yhg-de7119653_chrY_hg38full-woD_srt.bam':2.4555475863420606E-4):1.0000005000315726E-6):0.005724823759791775,dummy1):0.001, dummy2)0.01 ,dummy3 ):dummy4)dummy5:0.2)dummy6:0.3)dummy7);" ;

int main (int argc, char *argv[])
{
  FILE *f ;
  --argc ; ++argv ;		// swallow executable name
  if (!argc)
    { if (!(f = fopen ("ntest.in","w"))) die ("can't write ntest.in\n") ;
      fprintf (f, "%s", testString) ;
      fclose (f) ;
      if (!(f = fopen ("ntest.in","r"))) die ("can't read ntest.in\n") ;
    }
  else
    if (!(f = fopen (*argv, "r"))) die ("failed to open test newick file %s", *argv) ;
  TreeNode *n = readBinaryNewickTree (f) ;
  fclose (f) ;
  writeBinaryNewickTree (stdout, n) ;
  treeDestroy (n) ;
}

#endif

/*********** end of file *************/
