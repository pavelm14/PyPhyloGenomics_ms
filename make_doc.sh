pandoc -s -S -N -f markdown -V geometry:margin=1in  MS.md --bibliography=refs.bib  --csl=sysbio.csl  -o MS.pdf
