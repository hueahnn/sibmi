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
- O2 has been down :( sbatch still pending but see if the call works, next step is run with more files. To sample a random cluster of files from a directory: 
`ls <file directory> | shuf | head`

# Jun 20, 2025

- More on the background of the project and the so what:
  - Want to see the how prevalent heterodimerization due to insertion sequences are amongst all heterodimers
  - AMR genes historically were more prevalent on the chromosome but since the introduction of antibiotics have been moved to plasmids. 
- Breakdown of project in steps:
  1. Set up
  2. Annotations --> download all meta data and use pandas to filter through for only the complete and circular sequences.
  3. Find insertion sequences and ORIs
  4. Of the heterodimers, filter out ones with the specific topology of IS ORI IS ORI
- To do list for today: 
  - Get Snakefile to work for large datasets (having each bakta call submit a sbatch job to the cluster)
  - Set up pandas and Jupyter and filter through the metadata (and download the metadata onto the cluster)
    - Instead of using the OS package, use pandas to extract the names of the files
- Snakemake command: `snakemake --profile slurmprofile --rerun-incomplete --use-conda --executor slurm`
- For exporting a conda env to a JupyterLab kernel: `python -m ipykernel install --user --name=firstEnv`

# Jun 23, 2025

- Potential issues + cause of issues:
  - User was renamed at some point in time and that is screwing up the admin rights and read/write permissions
  - Also wsl may potentially be corrupted but uninstalling wsl would remove all the packages and other data associated with it
- Task for the week while O2 is down: create a workflow for the downstream analysis of the annotated sequences
  - Identify ORIs and IS using built software
    - OriV-Finder, ISEScan, ABRicate
  - Figure out a way to index the features so that only ones that have the correct topology of ORI IS ORI IS can be filtered out
  - Test on the 3 sequences pulled from literature that are confirmed to have the specific sequence we are looking for.
  - Make this into a snakefile
- Test batch 3 seqs from the papers
  - Done finding ORIs and ISes
  - How am I actually going to implement finding the correct pattern?
    - Key features: same IS seen 2+ times, ORI sandwiched between any of the same repeating IS

# Jun 24, 2025

- Abricate → used to find ORIs (for now)
  - Example usage: abricate --db plasmidfinder fusion-seqs/*.fasta > ORIs.tsv
- Approach to finding the heterodimers of our interest
  - Determine the number of class changes (don’t want any that have 2 or less)
- Steps:
  - First compile the results into one big .tsv file and do some filtering with pandas
    - Ex: filter out any IS that only show up once
  - Since all the positions exist, can determine the order of the features and look how many “class changes” occurs for each IS
  - Positive controls: the 3 plasmids from literature, neg control: the monomers from literature
- Snakefile to do:
  - Run abricate on all files to find ORIs
  - Run ISEScan on all files to find IS and compile all results into one .tsv file
  - Run heterodimer-finder.py script
- Application of project to antibacterial resistance
  - Can run plasmids through AMRfinder and can present a summary finding of there are this many plasmids that gained AMR genes through plasmid dimerization via IS

# Jun 25, 2025

- Finished heterodimer-finder.py
- Working on Snakefile
- ISEScan not working for some reason → troubleshoot
- Problem:
  - 2 of 3 positive controls successfully identified, however one failed because none of the IS identified in the paper were identified by ISEScan
  - Fix: BLAST search using ISEFinder database


# Jun 26, 2025

- Generalize heterodimer-finder.py to 2 scenarios:
  - Aggregated file with information about multiple sequences
  - For a single plasmid
- ISEScan not working… made a new conda env with only ISEScan and this seems to be working
- Snakefile fixed and working
  - 2 ways of finding IS: ISEScan or BLAST search using ISEFinder database
- Problem: one of the negative controls is being picked up as a sequence of interest
  - Not a problem → not uncommon to have multiple origins and could be sandwiched between ins seqs by chance
  - Also, we don’t know for certain this is a negative control. 
- Co-occurrence of ORIs → what’s rare what’s not
- 2-step BLAST: one of each sequence against the ISFinder database, then the result of each of that (the ones with multiple occurrences) extract the sequence and pairwise BLAST (because there may be subtypes)

To do:
- Annotations with bakta → o2
  - Fix Snakefile and start running on a larger set of seqs (maybe 100 or so)
- Pairwise BLAST → find subtypes for each IS
  - Output tsv file from initial BLAST has no headers… this is annoying → guide
  - Make tweaks to heterodimer-finder.py to support BLAST output
- Analyze co-occurrence frequency of identified ORIs 
- Install OriVFinder on o2
  - Change Snakefile accordingly

# Jun 27, 2025

Priority List:
- Submit annotation job with 100+ seqs
  - Create a script to parse through all seqs and filter only for complete and circular seqs. Also tsv only contains seqIDs, need to somehow obtain the associated fasta file.
  - Can obtain fasta file by appending .fasta to the end of the seqID and downloading the databases
  - Combine all the files from the databases into one directory and make sure they all have the same extension
  - Tweak Snakefile to read in seqs from txt file and test batch on SLURM with 100+ seqs
- Install OriVFinder on o2

To do:
- Testing 15 files right now, if successful launch job with 100 files and let it run over the weekend. 
# Jun 28, 2025
- 1000 worked successfully. Running with all ~75,000 however am missing some input seqs:
  - mMGEs_DRR003613_141_31316_r1256~5754148
  - mMGEs_DRR003613_141_31340_r1256~1018561  
  - mMGEs_DRR003614_141_36435_r1256~4396190
  - mMGEs_DRR003616_119_25658_r1256~3570481
  - mMGEs_DRR003616_119_25665_r1256~3570484

