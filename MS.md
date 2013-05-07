% PyPhyloGenomics: toolkit and protocol for developing phylogenetic markers in novel species for Next Generation Sequence data
% Carlos Peña^1^; Victor Solis^2^; Pável Matos^3^; Chris Wheat^4^
% 2013-05-01

^1^ Laboratory of Genetics, Department of Biology, University of Turku, Turku, Finland. Email: <mycalesis@gmail.com>

^2^

^3^

^4^

# Introduction
Next Generation Sequencing (NGS) is considered a quantum leap in improvement
in techniques for DNA sequencing.
The sequencing output of NGS technology is around 30 gigabases of DNA in one
single run [@reis2009] while the traditional Sanger method [@sanger1977] allows
sequencing only \~1,000 bp per specimen in the old capillary-based technology.
This higher yield is achieved by using massive parallel sequencing of PCR
products based on DNA synthesis using micron-scale beads on planar substrates
(a microchip) [@shendure2008].
As a result, millions of copies of sequences (reads) are produced from the DNA
templates. 
One application of NGS is targeted sequencing of numerous loci of interest
[@ekblom2010] in one run, which is quicker and cheaper than
using the traditional Sanger method.

Research in phylogenomics, using many more genes from genome sequencing
[@wahlberg2008], would be accelerated by using NGS due to the ease to obtaining
DNA data at massive scale. It would be very easy  to sequence many more than 
the 12 to 19 loci that so far have been used in phylogenomic studies 
[@wahlberg2008; @regier2013].

However, researchers have been relying on the Sanger method for sequencing a
handful of genes to be used in phylogenetic inference in several Lepidoptera
groups [@matos2013; @regier2013; @pena2011].

One issue to develop is a way to obtain candidate genes suitable for
phylogenetic inference, i.e. orthologous, single copy genes, lack of introns, etc.

* why orthologous
* why single copy
* why no introns
* why separated by xxxx distance

Indeed some studies have used NSG techniques to study phylogenetics at the 
genomic level using miRNAs for higher level phylogeny in Arthropoda
[@campbell2011]. miRNAs are nonprotein coding RNAs of small length involved in
DNA transcription and gene regulation. Using miRNAs for phylogenetics has the
drawback that these molecules are not easy to sequence from genomic DNA as miRNAs
are processed in the cell and shortened to \~ 22 base pair sequences
[@wienholds2005].

@regier2013 obtained their sequences from mRNA by performing reverse transcription
and PCR amplification [@regier2007]. mRNAs are molecules transcribed from
genomic DNA that have had introns spliced and exons joined. Therefore, attempting
to sequence these genes from genomic DNA for other species will be troublesome
due to the likely appearance of introns. Intron sequences can be of various
lengths across taxa and would prove difficult to assess homology for phylogenetic
studies.

@wahlberg2008 obtained candidate genes for phylogenomics by identifying single
copy and orthologus genes of *Bombyx mori* from EST libraries. EST sequences
were searched in the *Bombyx mori* genome in order to identify suitable exons.
These exon sequences were compared against EST libraries of related Lepidoptera
in order to obtain homologous sequences for primer design. Thus, this method
depends on the availability of EST sequences which are single reads of cDNA that 
might contain numerous errors and are prone to artefacts [@parkinson2002].

According to @wahlberg2008, it is easier to employ  genomic DNA for phylogenetic
practice due to several reasons: (i) genomic DNA does not degrade so quickly as RNA;
(ii) it is simpler to preserve in the field; (iii) it can be sequenced even from
dry material (for example museum specimens); and (iv) it is the most commonly used
DNA in molecular systematics.

One strategy to fulfill this goal could be comparing genomic sequences of  model
species and extract suitable genes that can be sequenced in novel species from its
genomic DNA. 

In this paper, we describe a protocol for finding genes in genomic DNA suitable for
phylogenomic studies. We have developed bioinformatic tools to help harvesting genes
from  genomes available in public databases. We have also optimized a wet lab protocol
for sequencing the found genes using NGS technology. Additionally, we  have developed
bioinformatic tools for analysis of the raw data from NGS. Our software has been
developed to filter the output reads from NGS and assemble the sequences for each
specimen and their sequenced genes so that datasets can be assembled for analysis in 
common software for phylogenetic inference.


# Methods

## Finding candidate genes from *Bombyx mori*

We decided to use the *Bombyx mori* genome as starting point (although any genome can
be used) to obtain candidate genes suitable for sequencing across novel species.
As explaine in the introduction, genes to be used in phylogenetic inference have to
fulfill the following requirements: (i) Our genes should be orthologs; (ii) Our genes
should be single-copy genes; (iii) their sequence need to be around 251 DNA base pairs
in length for easy sequencing in our Next Generation Sequencer available, IonTorrent_.

The OrthoDB database <ftp://cegg.unige.ch/OrthoDB6/> has a catalog of orthologous
protein-coding genes for vertebrates, arthropods and other living groups.

We parsed this list with the module OrthoDB from our package PyPhyloGenomics and 
obtained a list of single-copy, orthologous gene IDs for *Bombyx mori* (12 167 genes
in total).

By using the module BLAST of PyPhyloGenomics, we were able to obtain the gene sequence
for our list from the CDS sequences for *Bombyx mori*. And this list of sequences
were blasted against the *Bombyx mori* genome in order to make sure that our sequences
are free of introns.

We also filtered these genes to be at least 300bp in length, are separated by at least 
810kb so that they can be considered independent evolutionary entities.

At the end we obtained 575 exons in frame.

Then we executed exon validation by searching these exons in other genomes of Lepidoptera
species such as, *Danaus*, *Heliconius* and *Manduca*,

This search is automated by using functions in our module BLAST that take as input 
the list of genes from *Bombyx mori*, and the file with genomic sequences for the target
species.

PyPhyloGenomics also has tools to create fasta files for those exons that were found in
the other species' genomes. It also automates the alignment by using MUSCLE.

PyPhyloGenomics contains functions to automatically design degenerate primers from the
homologous sequences by delivering the sequences to primer4clades and receiving the 
designed primers.

It is recomended that both alignment and designed primers to be analyzed carefully 
to make sure that the are no problems.

We have to consider the requirements for the IonTorrent platform 2 in order to arrive
to our target 250bp gene length (Table 1).


For IonTorrent <http://www.iontorrent.com/> Platform 2, the maximum length that can be sequenced is from 280bp to 320bp in total. Thus, ``320 - 119 = 201`` is the maximum internal gene region (region within degenerate primers).

Therefore, for the new set of primers, being designed for Platform2, we have a maximum amplicon size of ``201 + 25*2 = 251bp``. 

**Table 1.** Adaptors and primers needed for sequencing in the NGS Iontorrent platform 2. 
The maximum length of sequenced amplicon is \~ 201 bp after discarding primer regions.

  Primer                Length (bp)
  --------------------  ------------
  Adapter A             30
  5' Index              8
  5' Degenerate Primer  25
  Exon                  **???**
  3' Degenerate Primer  25
  3' Index              8
  Adapter P             23
  --------------------  ------------



## Other software
CEPiNS [@hasan2013] is a software pipeline that uses predicted gene sequences
from both model and novel species to predict and identify exons suitable for
sequencing useful for phylogenetic inference.







## Action items
### Single Copy Genes
From OrthoDB and get a list of single-copy genes for *Bombyx mori*
Get the exon sequences from the CDS file for *Bombyx*
Use these CDS sequences to extract the corresponding sequence in the *Bombyx* genome avoiding gaps
so that we will only work with genes or gene fragments that do not include introns (which is good 
if you want to do phylogenetics).

## Exon Structure
Comparison against other butterfly and moth genomes 
Exon validation against genomes of *Manduca*, *Danaus* and *Heliconius*
 
## Exon Aligment and Primer Design
Need to get at least 4 sequences for primer desing using primers4clades [@contreras2009].

## Wet Lab
### Multiplex PCR
### Sample preparation for Next Generation Sequencing in IonTorrent

## Next Generation output analysis
### find barcodes
### find primers
### quality control of reads
### de novo assembly using velvet
### alignment and storage in VoSeq?
So that we can manage the high number of sequences and create datasets for phylogenetic analysis
very easily.

## Comparison with other methods
@regier2009, @regier2008, @regier2013 use Reverse Transcription PCR from mRNAs to avoid sequencing introns, although the corresponing genomic DNA sequences are  likely to include introns. Therefore if one use their genes, it is not recommended to do "direct gene amplification" [@regier2007].

# References
Targeted sequencing [@godden2012]




