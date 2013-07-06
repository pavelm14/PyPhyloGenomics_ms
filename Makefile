MS.pdf: MS.md refs.bib
	pandoc -s -S -N -f markdown -V geometry:margin=1.5in MS.md --bibliography=refs.bib --csl=style/molbiolevol.csl -o MS.pdf
