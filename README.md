# Finding TOTOs--SIBMI 2025 project with the Baym Lab.

## Contents
Main directory (`/n/data1/hms/dbmi/baym/hue/sibmi`):
1. analysis
  - `plots.ipynb`: reads in final dataframe to generate plots and analyze contents
  - `Snakefile`: contains rule amrfinder
2. bakta-annotations: contains all bakta annotations
3. conda-envs
4. fasta-files: contains all fasta files from the PlasmidScope DB
5. heterodimer
  - `blast.py`: parse ins seq BLAST output table + run pairwise BLASTs
  - `heterodimer-finder.py`: find TOTOs + append additional data to generate a final df
  - `heterodimer-plasmids-all.txt`: txt file containing the plasmid IDs of TOTOs (type 1 ORIs only)
  - `how_many.py`: figuring out how many and which files are empty
  - `origins.py`: all things ORI related-- find mORIs, stitch together orivfinder outputs, assign ORI, functions for co-occurrence 
  - `Snakefile`: finds ORIs and IS
- `Snakefile`: for bakta annotations