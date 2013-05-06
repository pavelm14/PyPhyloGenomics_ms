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
(ii) it is simpler to preserve in the field; (iii) it can be sequenced even
dry material (for example museum specimens); and (iv) it is the most commonly used
DNA in molecular systematics.

Thus, we have developed a protocol for finding genes in genomicd  DNA suitable for
phylogenomic studies. We have developped bioinformatic tools to help harvesting genes
from  genomes available in public databases. We have also optimized a wet lab protocol
for sequencing those found genes using NGS technology and more bioinformatic tools
for analyisis of the raw data from NGS. Our software has been developed to assemble
the wanted sequences from the reads of the NGS machine so that datasets ready for
anlysis in phylogenetic inference can be created.


## Other software
CEPiNS [@hasan2013] is a software pipeline that uses predicted gene sequences
from both model and novel species to predict and identify exons suitable for
sequencing useful for phylogenetic inference.





## Exon models 	
### Finding candidate genes from *Bombyx mori*

We need to obtain candidate genes to be used in phylogenetic inference that have to fulfill the following requirements:

* Our genes should be orthologs.
* Our genes should be single-copy genes.
* Their sequence need to be around 251 DNA base pairs in length.

We will assume that our Next Generation Sequencer available is the IonTorrent_.

We have to consider the IonTorrent_ platform requirements to arrive to our target 250bp gene length:

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

For IonTorrent <http://www.iontorrent.com/> Platform 2, the maximum length that can be sequenced is from 280bp to 320bp in total. Thus, ``320 - 119 = 201`` is the maximum internal gene region (region within degenerate primers).

Therefore, for the new set of primers, being designed for Platform2, we have a maximum amplicon size of ``201 + 25*2 = 251bp``. 

The OrthoDB <ftp://cegg.unige.ch/OrthoDB6/> database has a catalog of orthologous protein-coding genes for vertebrates, arthropods and other living groups.


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




