% PyPhyloGenomics: toolkit and protocol for developing phylogenetic markers in novel species for Next Generation Sequence data
% Carlos Peña^1^; Victor Solis^2^; Pável Matos^3^; Chris Wheat^4^
% 2013-05-08

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

Research in phylogenomics would be accelerated by using NGS due to the ease
to obtaine DNA data at massive scale. It would be very easy  to sequence many
more than the 12 to 19 loci that so far have been used in phylogenomic studies 
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
copy and orthologus genes of *Bombyx mori* from EST libraries. They searched for
EST sequences in the *Bombyx mori* genome in order to identify suitable exons.
These exon sequences were compared against EST libraries of related Lepidoptera
species
in order to obtain homologous sequences for primer design. Thus, this method
depends on the availability of EST sequences which are single reads of cDNA that 
might contain numerous errors and are prone to artefacts [@parkinson2002].

According to @wahlberg2008, it is easier to employ  genomic DNA for phylogenetic
practice due to several reasons: (i) genomic DNA does not degrade so quickly as RNA;
(ii) it is simpler to preserve in the field; (iii) it can be sequenced even from
dry material (for example museum specimens); and (iv) it is the most commonly used
DNA in molecular systematics.

Thus, it is needed a way to find candidate genes that can be easily sequenced
from genomic DNA across several lineages.
One strategy to fulfill this goal could be comparing genomic sequences of  model
species and extract suitable genes that can be sequenced in novel species from 
simple extractions of genomic DNA. 

In this paper, we describe a complete protocol for finding genes in genomic DNA
suitable for phylogenomic studies, sequencing of those genes using NGS techonology,
and data analysis of the NGS output.
We have created the software package PyPhyloGenomics that includes bioinformatic
tools required for automated gene finding, primer design and NGS data analysis
Our software has been
developed to filter the output reads from NGS and assemble the sequences for each
specimen and their sequenced genes so that datasets can be assembled for analysis in 
common software for phylogenetic inference.


# Methods

## Finding candidate genes from *Bombyx mori*

We are interested in studing the phylogenetic relationships of lineages in the 
Lepidoptera. Hence, we decided to use the *Bombyx mori* genome as starting point 
(although any genome can be used) to obtain candidate genes suitable for sequencing
across novel species.
As explained in the introduction, genes to be used in phylogenetic inference have to
fulfill the following requirements: (i) the genes should be orthologs; (ii) the genes
should be single-copy genes; (iii) their sequence need to be around 251 DNA base pairs
in length for easy sequencing in our in-house Next Generation Sequencer, an IonTorrent
Ion PGM Sequencer from Life Technologies (<http://www.iontorrent.com/>).

The OrthoDB database <ftp://cegg.unige.ch/OrthoDB6/> has a catalog of orthologous
protein-coding genes for vertebrates, arthropods and other living groups.
We parsed this list with the module OrthoDB from our package PyPhyloGenomics and 
obtained a list of single-copy, orthologous gene IDs for *Bombyx mori* (12 167 genes
in total).

A function in our module BLAST extracted the sequences for those genes from
the *Bombyx mori* CDS sequences (available at <http://silkdb.org>). We BLASTed
the sequences against the *Bombyx mori* genome and discarded those containing
introns.
We kept genes with sequences longer than 300bp in length, and separated by
at least 810kb so that they can be considered independent evolutionary entities and
obtained 575 exons.

We validated those exons by automated search  of these exons in other genomes of
species in Lepidoptera, such as, *Danaus*, *Heliconius* and *Manduca*.
This search is automated by using functions in our module BLAST that take as input 
the list of genes from *Bombyx mori*, and the file with genomic sequences for the
target species.
During validation, PyPhyloGenomics creates FASTA format files by appending
matching sequences from the tested genomes. It also automates the alignment of
sequences by using the software MUSCLE.

PyPhyloGenomics contains functions to automatically design degenerate primers
from the homologous sequences by delivering the sequences to primer4clades
and receiving the designed primers. primer4clades is a web service based on the
CODEHOP strategy for primer design [@contreras2009].

It is recomended that both alignment and designed primers to be analyzed carefully 
to make sure that the are no problems. After this step, one can have around XXX
genes ready to be sequenced across novel species in Lepidoptera for many species
if NGS techniques are used.

### Sample preparation for Next Generation Sequencing in IonTorrent
We followed the library preparation protocol for NGS by @meyer2010 that can be 
adapted for the IonTorrent technology. This method consist in attaching and index
(or barcode) to the sequences of each specimen for later identification. Therefore,
it is possible to separate reads from the NGS data using these indexes.

We sequenced several individuals of a wide range of species in the Lepidoptera. We
also sequenced a specimen XX_Code_XX of *Bombyx mori* as control for our lab and
NGS assembly data protocols.

The IonTorrent platform 2 can sequence from 280 to 320bp in total. The Iontorrent
adapter, index and primer sequences leave around 201bp as the maximum internal
gene region sequenced (region within degenerate primers) (Table 1).

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

## Other software
CEPiNS [@hasan2013] is a software pipeline that uses predicted gene sequences
from both model and novel species to predict and identify exons suitable for
sequencing useful for phylogenetic inference.

# References
Targeted sequencing [@godden2012]




