# Jun 9, 2025

Rundown of project with Kepler on the M2
- Plasmids with insertion sequences are cool in that they can combine and form fusion plasmids (dimers) and this is a way for antibiotic resistance genes and other mutations to spread. Essentially they form this one big plasmid.
- Can identify these plasmids because it will have 2 ORIs. Specifically IS ORI IS ORI
- Find plasmids with this pattern (2+ repeating) 


# Jun 11, 2025

Instructions from Kepler:
1. Set up cluster
   1. Get O2 account
   2. Install conda/mamba
   3. Install jupyterlab
2. Download the PlasmidScope database → done by Arya, on o2
3. Homogeneously annotate the PlasmidScope database
   1. Use snakemake with bakta
   2. If this is too much only run it on circular complete plasmids

Getting set up with Arya:
- Most important: getting set up with o2
  - Note: o2 will be down 6/23-26… use git to access code and work locally
  - github repo: https://github.com/hueahnn/sibmi 
  - Go through a git repo tutorial 
- Snakemake: essentially a directed acyclic graph. Use defines inputs that match to an output and a function maps each input to an output. This “function” is user defined script
  - Documentation: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html 
  - Example: https://github.com/karel-brinda/MiniPhy/ 
  - Read through snakemake documentation
- [PlasmidScope](https://plasmid.deepomics.org/download#gbk) Database
  - About: https://plasmid.deepomics.org/tutorial 
- bakta annotation software: https://github.com/oschwengers/bakta
  - bakta DB path on o2: /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db
  - example usage: bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db ../../decompressed_genomes/allthebacteriav2/stdin.part_002/nisin-susceptibility-associated-two-component-system-response-regulator-NsaR_GCF_900635095.1_33763_A01_WP_000697886.1-len666/SAMN18961669.fa --output SAMN18961669_bakta --force --skip-plot
- More down the line, but familiarize with pandas if using python: https://pandas.pydata.org/ 

Tips + words of advice:
- Know why I’m doing what I’m doing and have a clear narrative and story behind my project for the end of the summer presentation
- For the snakemake pipeline I will be building, do it in batches and scale up → good practice
- As a bioinformatician, it is easy to get lost in the sauce 🥫make sure to understand the biology behind what I’m doing and ask lots of questions!


# Jun 12, 2025

Action items for today: 
- Git clone and set up repo
- Make snakemake pipeline

Snakemake pipeline: goal – to annotate all of the plasmids in the database using bakta

Familiarize with:
- Git version control and documentation
- Snakemake
- Clusters and shared compute systems

Tomorrow:
- Bakta trial
- Snakemake pipeline


# Jun 13, 2025

- Snakefile time → can just use vim
- Example bakta run: 
`bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db ./PLSDB/fasta/PLSDB_OX335189.1.fasta --output ./bakta --force`
- Git set up
  - Use issues to document → kind of like a to-do list/diary
  - Use branches to keep neat
    - Branches open: master, snakemake
    - git checkout -b <name>
- Requesting an interactive o2 job: 
`srun --pty -p interactive -t 12:00:00 --mem=16G bash`
- Extracting + unzipping files: 
`wget <url>`
`tar - xvzf <file>`


# Jun 16, 2025

- Snakefile works for one file → now generalize for multiple files and scale up
  - Points of confusion: 
    - If calling on multiple files, how to do this because with one file I specify the name of the output file when I call snakemake in the command line
    - rule all structure: check each folder to see if all output files are present and if so output to a .txt file “complete”
- Works for one file, generalized to 2 then 8
  - When running 8 got a memory error → either increase the memory (>16G) or decrease the core (was running on 4, decreased to 2 which was much slower but ran through the end)
  - 8 complete → but took a while, generalize to larger batches by having each call submit a new sbatch job
- Arya sent his code for submitting cluster jobs within a snakemake call

# Jun 17, 2025

- Implementing cluster job submission feature, to call:
`snakemake --use-conda --profile slurmprofile -p --keep-going --rerun-incomplete --executor slurm`



