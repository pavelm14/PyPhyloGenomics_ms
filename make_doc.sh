pandoc -s -S -N --template header.latex -f markdown -V geometry:margin=1in  MS.md --bibliography=refs.bib  --csl=style/molbiolevol.csl  -o MS.pdf
