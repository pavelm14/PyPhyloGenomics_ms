# This is the procedure to reproduce the analyses in the PyPhyloGenomics manuscript
#

## Get data files
# Orthodb file
data/OrthoDB6_Arthropoda_tabtext.csv: 
	wget ftp://cegg.unige.ch/OrthoDB6/OrthoDB6_Arthropoda_tabtext.gz
	mkdir -p data
	mv OrthoDB6_Arthropoda_tabtext.gz data/.
	gunzip data/OrthoDB6_Arthropoda_tabtext.gz
	mv data/OrthoDB6_Arthropoda_tabtext data/OrthoDB6_Arthropoda_tabtext.csv

# Bombyx mori CDS
data/silkcds.fa:
	wget ftp://silkdb.org/pub/current/Gene/Glean_genes/silkworm_glean_cds.fa.tar.gz
	mv silkworm_glean_cds.fa.tar.gz data/.
	tar -zxvf data/silkworm_glean_cds.fa.tar.gz -C data/

# Bombyx mori genome 
data/silkworm_genome_v2.0.fa: 
	wget ftp://silkdb.org/pub/current/Genome/silkworm_genome_v2.0.fa.tar.gz
	mv silkworm_genome_v2.0.fa.tar.gz data/.
	tar -zxvf data/silkworm_genome_v2.0.fa.tar.gz -C data/

gene_search: data/OrthoDB6_Arthropoda_tabtext.csv data/silkcds.fa data/silkworm_genome_v2.0.fa
	echo "Done"

