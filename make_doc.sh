pandoc -s -S -N -f markdown -V geometry:margin=2in  MS.md --bibliography=refs.bib  --csl=sysbio.csl  -o MS.pdf
