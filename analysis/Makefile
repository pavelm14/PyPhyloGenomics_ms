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

gene_search: data/OrthoDB6_Arthropoda_tabtext.csv
	echo "Done"

