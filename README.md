# phynder

Efficient likelihood calculations to place samples into a phylogenetic tree.  In particular, phynder was originally designed for placing male ancient DNA samples into the Y chromosome phylogeny, when there are arbitrarily high rates of missing data.  It has also been used for assigning ancient samples to whole mitochondrial genome phylogenies.

### Synopsis

```
phynder -q query.vcf -p 0.01 -o query.phy tree.nwk tree.vcf.qz
```
where query.vcf contains one sample called **test**. This gives the
following output in query.phy

```
1 3 smp 1 1
! 4 7 phynder 3 1.0 62 phynder -q query.vcf -p 0.01 -o query.phy tree.nwk tree.vcf.gz 19 2020-05-14_13:21:02
.
~ D I 1 6 STRING                sample name
~ D B 3 3 INT 4 REAL 4 REAL     best branch, posterior, score
~ D C 1 3 INT                   clade
~ D N 1 3 INT                   number of sites used to assign branch
~ D S 3 3 INT 4 REAL 4 REAL     suboptimal branch, posterior, score
.
. read tree with 4027 nodes and 2014 leaves
. read 121335 sites for 2014 samples in tree
.   ignored 631 multiple-allelele and 0 non-SNP sites
. built scores for 121323 gt patterns
. using 120908 good sites
.   427 sites rejected because likelihood below threshold -10.0
.   121335 sites had missing genotypes, mean 31.5
. read 7221 sites for 1 samples in query
. 
I 4 test
N 7154
B 424 0.478123 -54.000068
S 413 0.020019 -57.173241
S 414 0.034168 -56.638637
S 419 0.007578 -58.144716
S 423 0.449961 -54.060776
C 409
```
This is a ONEcode file, which is a general data self-documenting data
exchange format in which each file is composed of a header made of lines
starting with symbols that are not letters, followed by a body
containing data made of lines starting with letters.  See the [vgp-tools website](https://github.com/VGP/vgp-tools) for
further details on the general ONEcode format and associated tools.

The first line of the header
identifies this as ONEcode, defines the file type (smp for "sample")
and format version.  The next line provides the provenance of the
file. The next block gives the schema, defining each line type that
will be used in the body of the file.  Then come a set of comment
lines that in this case describe some key facts of the data
processing.

The body is composed of lines that start with alphabetic characters,
following the schema. In our case I lines give the sample identifiers
and define the objects,
B lines give the best branch assignment with its posterior probability and
likelihood score, and S lines give suboptimal branch assignments with their
posteriors and scores.  The output lists suboptimal assignments with
posterior probability higher than the *p* parameter given on the
command line, in this case 0.01. The C line indicates the the most
specific branch in the tree that contains at least (1-*p*) of the total
posterior probability.

```
phynder  -B -o branches.snp tree.nwk tree.vcf.gz
```
Here the file `branches.snp` will be created, listing for each of the
variant sites in `tree.vcf.gz` the most likely branch in `tree.nwk` on
which the relevant mutation occurred, plus some other information described below.

### Installation

Phynder requires htslib, the library that reads and writes SAM/BAM/CRAM and VCF/BCF files.
To obtain phynder:

```
git clone https://github.com/samtools/htslib.git
cd htslib
make
cd ..
git clone https://github.com/richarddurbin/phynder.git
cd phynder
make
make install
```

NB `make install` copies the binary into your ~/bin directory.  If you want to install it elsewhere, either
copy the executable there yourself or change the makefile to copy to where you want.

### Commands and options

To see the command line options run `phynder` without arguments.  You should get a help message with
the available options/commands.

The general form of a phynder command line is a series of options, followed by two required file names. 

```
  phynder [options] <newick tree> <vcf for tree>
```
Further information on options is given below.  They can come in any
order, because they are all parsed before any computation is carried out.

The tree must be in standard Newick format, and should have branch
lengths corresponding to times, ideally built with a maximum
likelihood method using the data in the VCF file. It can either be
*rooted* with two branches deriving from the top level node, or
*unrooted* with three branches deriving from the top level node.

The sites that will be used in the analysis must have genotypes given
for the samples leaves in the tree in a VCF file, which can in
practice be an uncompressed or gzip compressed .vcf file, or a binary
.bcf file as made for example with bcftools.  Only biallelic single
nucleotide variants are used.  The user does not need to assign the
ancestral state for variants.  Analysis will be faster if there is no
missing data, but some level of missingness is tolerated.  High
missingness will increase compute time and degrade performance.

```
  -o <filename>            output filename - this will be a ONEcode file
  -b                       output in binary
```


```
  -ts <transition rate>    [1.3300]
  -tv <transversion rate>  [0.3300]
```
These set the mutation rates for transitions and transversions.  The default values are appropriate for the human Y chromosome if the branch lengths in the tree were built with an expected mutation
rate of one mutation per site.

```
  -T <thresh>   site log likelihood threshold [-10.0]
```
This sets a threshold for the log likelihood a site must achieve to be used in subsequent analysis.  The likelihood of a site can be too low if there were many independent mutations at the site, but
more likely because the site is prone to error in genotype calling. The units are natural log
units. The true threshold is in fact this value plus the log of the probability of a mutation on the
shortest branch, since we need to allow for a single mutation on any branch.

Setting `T 0` will disable the theshold, using all sites.  This may in particular be of interest in conjunction with `-B` so as to collect information on all sites.  Otherwise the default value is in most cases sensible.

```
  -q <query vcf>           a vcf file of samples that the user wants to place into the tree
  -C <calc_mode>           [0] calculation mode
                           calc_mode 0: LL both ends of edge match
```
The -q option that places additional samples into the tree.  The query vcf can contain multiple
samples.  For each sample it will output a line with the best branch position, its posterior
probability, the number of sites used in that calculation, the total score and the score per site.  A score per site close to 0 indicates good overall fit of the sample data to the tree.  Something is wrong it this is below -0.1.

The `-C` option is present because the calculation used by phynder is, for reasons of efficiency, an approximation.  The sample will in fact branch off the middle of some edge in the reference tree and mutations actually on that edge are informative about where on the edge to place the junction. However, the default (calc\_mode 0) behaviour of phynder fails to correctly use this information, though it should correctly assign the edge. (Formally, it calculates the posterior conditional on the whole incident edge matching the query.)  There are various other possible ways to do this, with different degrees of computational efficiency, which could be implemented as alternative options (see "To do" section below). 

```
  -p <posterior threshold>  give suboptimal branches and clade [0.]
```
This option must precede the `-q` command, and if present and given a non-zero value generates out two extra types of output line per query.

First, it prints out any suboptimal branches with posterior probability higher than the given value *p*.  

Second, it finds the lowest branch in the tree for which the sum of posterior probabilities for the whole clade below that branch (including the branch in question) is greater than (1 - *p*).  This provides a conservative assignment when placement is uncertain.  

When data are sparse, there may be a whole subregion of the tree for which no genotype information is informative, and in this case longer branches get higher posteriors just because they are longer.  In this clade placement is probably more appropriate than just picking the best fit branch.

```
  -B <branch file name>    output branch positions of tree variants
```
This writes out a file assigning each site in the VCF file (that passes the threshold) to a branch.
The output columns are: CHRomosome, POSition, REF allele, ALT allele, isANCestral, BEST edge,
BESTSCORE, NEXT edge, NEXTSCORE, LOGLIKE, where ANC is 1 if the ALT allele is ancestral, the
BESTSCORE is the log likelihood of the data above the edge being ancestral, and the data below the
edge being derived, and LOGLIKE is the log likelihood of all the data at the site, including the
probability of a mutation on the best edge.  This means that LOGLIKE is often lower than BESTSCORE.

Like `-q`, the `-B` option must follow `-t` and `-v`.

```
  -h       print the help message
  -V       verbose - print extra info
```
These are simple utility commands.  `-V` gives more information on missing data in the reference vcf, timing and memory use information of various steps, and some other information.

```
  -U   make tree ultrametric - all leaves equidistant from root
```
In principle trees generated from modern samples under a uniform clock should be ultrametric, so there is some rationale to enforce this, but it is not clear that the benefits outweigh the loss of fit due to changing the branch lengths from the maximum likelihood values.  In most cases the differences in outcomes should be very small.  However this option will be required in order to use a coalescent prior, when that is implemented.

If used this would need to come before the `-B` or `-q` commands.

### To do
   
* Add an option to use a coalescent prior.  This is not so hard, and will be computationally cheap. Currently a uniform prior is used over branches. Note that because of the likelihood calculation approximation (see above) this does not mean that the posteriors will be uniform given a lack of informative data.

* Add the ability to read genotype likelihoods (ideally) and/or posteriors from the query VCF file. This is the right way to handle uncertainty in the query genotype calls. Currently phynder assumes calls are Q40 (probability of error 1 in 10,000).  Although this will be incorrect for low coverage data, particularly at deamination sites, it is not very likely that errors will be sufficiently consistent with the tree to be misleading, so the consequences of this assumption will generally be benign with respect to the result.
	* I could alternatively, or as well, allow the user to set generic likelihoods based on the mutation class, e.g. setting C to T transitions as more likely errors.
* Try alternate approximations for the calculation mode.

Richard Durbin
May 2020
