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
  - mMGEs_DRR003622_141_101836_r1256~2120518


# Jun 30, 2025

- mMGE missing file problem is because the directory name doesn't match the input names
- bakta annotations are running --> will take a while to finish running
  - had to increase memory to 24G bc some of the jobs were failing
  - instead of submiting one bakta per job, find a way to submit multiple baktas per job so that the queue doesn't get clogged up like today (the higher beings of o2 are pissed off rn ☹️). 

TODO:
- Run initial screen using abricate and ISEScan to get a preliminary dataset
  - Fix Snakefile to deal with the scenario in which ISEScan doesn't find any IS and returns no output files (workaround is expanding the directory in rule all)
- Install OriVFinder and update Snakefile accordingly 
  - asked o2 ppl for help, awaiting response
- Script for pairwise BLAST and figuring out which IS are pairs/unique
  - collapse duplicates (take the highest bitscore for seqs with multiple hits at the same location). filter out high e-vals. run pairwise BLAST searches on unique IDs


# Jul 1, 2025

- Heterodimer Problems:
  - BLAST broken... look into this (.blast.tsv files are empty)
  - Figure out how to skip seqs without insertion sequences-- ISEScan returns no output files for these
- Lower priority item: parsing through initial BLAST output, create a script to run pairwise BLASTs

Worklog:
- Ran 1000 bakta jobs, running another 1000 rn that is stuck in the queue
- Fixed heterodimer-finder.py script, running on 100 seqs

TODO:
- Keep running bakta jobs + keep running heterodimer-finder.py
- Start working on BLAST script while jobs are running tomorrow.
- How to look for plasmid fusion events --> looking for plasmids comprised of 2 distinct plasmids 
  - Lineage info may already exist, look into this
  - Can also BLAST, depends on how many plasmids get filtered out
- Slack lab person in Norway who studies plasmids


# Jul 2, 2025

- Bakta annotations are halfway done, ran overnight and currently at 35064 
- Out of 1000 seqs passed through heterodimer-finder only 22 have the pattern we are looking for --> although this is using ISEScan

TODO:
- Prelim results from 1000 random seqs:
  - How many have multiple origins and chart of co-occurence of origins
- Harvard geneious license

- Turns out abricate did not do a very good job at finding the origins... wait for OriVFinder to be installed.
  - 22 seqs had the wanted feature, only 46 of the 1000 seqs had 2+ origins... should be around 1/3


# Jul 3, 2025

TODO:
- blast.py and ORI co-occurence table
  - for the ORI co-occurence table need 2 things rn, frequency chart of how often 2 ORIs occur together and a freq chart of how often an ORI appears.
later:
- dual boot ubuntu
- Use Docker to run OriVFinder locally


# Jul 7, 2025

- For blast.py, use the plasmid seq as the db and use each insertion sequence hit as the query. From this filter for percent identity > 90
- 9 failed bakta jobs (PILER-CR error, error code: -11): skipping CRISPR for these ones...
  1. mMGEs_k99_473275_r1256~1033030
  2. COMPASS_KY595967.1
  3. GenBank_CP082893.1
  4. COMPASS_NC_020243.1
  5. PLSDB_NZ_LR890507.1
  6. GenBank_KY000037.1
  7. GenBank_CP081352.1
  8. GenBank_CP066478.1
  9. GenBank_CP053913.1
- Prioritize working on Docker, need preliminary results. If Docker is too much of a hassle, use MOB finder
- Trying to figure out Docker and OriVFinder is taking years off my life... fully crashing out😵‍💫
  - To run OriVFinder: `docker run -it --rm orivfinder-ready bash`
  - To execute command: `python oriVfinder.py --fasta /app/data/input.fasta --output_dir /app/data/output`


# Jul 8, 2025:

- Can run OriVFinder on O2 with some tweaks --> already have bakta annotations so skip that step to save time
- Try running OriVFinder locally using Docker to see what the output looks like
- Integrate BLAST script to Snakefile and run it for all the plasmids (might want to comment out some of the rules)
- Start thinking about final presentation and creating a compelling story.
  - Take a look at the paper Kepler sent.
- Order of business: get Snakefile with BLAST up and running, OriVFinder, have presentation in the back of mind.


# Jul 9, 2025:

- Rewrite blast.py, create a fasta file of all the unique hits as the DB for the BLAST search and query each individual fasta file 
- Check out PlasAnn output and determine which one to use as the ORI finder
- Good test case of insertion seqs: DDBJ_AP017321.1
- Path to OriVFinder on O2: `/n/data1/hms/dbmi/baym/kepler/orivfinder/app`
Progress:
- Running blast.py on snakemake
- Decide on OriVFinder vs PlasAnn for finding ORIs --> once this is decided can fix heterodimer-finder.py and get preliminary data. yay!! (also permitting that there are no issues with the Snakefile)


# Jul 10, 2025:

- Work on OriVFinder for O2
- Good examples of pairwise BLAST outputs: 
  - heterodimer/heterodimer/DDBJ_AP019793.1.pairwise.blast.tsv
  - heterodimer/heterodimer/DDBJ_AP019544.1.pairwise.blast.tsv
- Circle plot of association for ORIs


# Jul 11, 2025:

- Order of business:
  - Run OriVFinder + BLAST ins seqs on all 75000 plasmids
  - Figure out how many of the 1000 that already have run have non-empty files (for BLAST check if there are any outputs, for ORIs check how many type 1, type 1 or 2, and if there are multiple ORIs)
  - Parse BLAST pairwise output (group ins seqs together)
  - Update heterodimer-finder.py with new output formats
  - Make presentation


# Jul 14, 2025:
- Corrupted bakta annotation: COMPASS_NC_002664.1
  - Removed from txt file for now
- Another one corrupted: RefSeq_NZ_CP093161.1, this one is also a linear seq...
- Problems:
  - original tsv has a mix of tab separated and comma separated values
  - BLAST is not hitting all transposases (532 of 1000 do not have transposases)
- Make slides, Sneha visiting tmr
- To do:
  - fix origins.py
  - fix heterodimer-finder.py


# Jul 15, 2025:
- fixed and ran heterodimer-finder.py --> out of 999 seqs, 224 have desired feature
- `O2_jobs_report --report` to view memory usage and such


# Jul 16, 2025:
- How many have multiple ORIs? How many actual plasmids?
  - 365 with multiple ORIs, 100 plasmids with desired feature, 224 hits across these 100
  - Upper bound, haven't filtered the transposons. As a sanity check later on, plot lengths vs abundance and make sure that identified fusions have a longer length
- Name proposal: TOTO
Downstream analysis directions:
- Determine cooccurence freq of ORIs to filter results for seqs that are actually meaningful (ORIs frequently coexisting is not interesting)
  - Pull out ORI seq for every plasmid with 2+ ORIs and cluster at thresholds (70,80,90%) using mmseq2
  - Build a heat map with every permutation of ORIs by searching for the existence of the pair in every plasmid
- Are there any particular ISs that show up often?
- How many of these are events in which AMR genes (AMRfinder) or mobilizable genes (MOBfinder) are present?
To do:
- clustering with mmseqs2 --> not sure how to interpret the output
  - first col is cluster (representative seq) and second col is all seqs
  - 90% --> 717 clusters, 70% --> 696 clusters
  - prioritize type 1 hits only and cluster by 30% identity (ORIs are well conserved so little variance is to be expected)
  - to determine the actual identity will prob need to BLAST against OriVFinder DB


# Jul 17, 2025:
- 545 ORIs of type 1, clustering @ 30% identity --> 389 clusters
- Commands to run mmseqs2:
  - first convert FASTA DB to MMseqs2 DB: `mmseqs createdb examples/DB.fasta DB`
  - `mmseqs cluster DB DB_clu tmp --min-seq-id num`
  - `mmseqs createtsv DB DB DB_clu DB_clu.tsv`
  - to determine number of unique clusters: `import pandas as pd` `column_names = ['cluster', 'origin']` `df = pd.read_csv("___.tsv", sep="\t", names=column_names)` `df['cluster'].nunique()`
- BLAST ORIs against OriVFinder DB, assign top hit to each seq
- Make some charts with length, presence/absence of MOB & AMR genes
Future steps:
- Generate a chart of coocurrence frequencies and generate a circle plot of association
- Calculate the frequency of heterodimerization events
- Plot lengths vs abundance
- Create bar graphs comparing AMR, MOB, # of IS/per kb for 2+ ORIs vs TOTOs
- Future project --> look for heterodimers via kmers, can use kmer density to determine dimers formed by insertion seqs


# Jul 18, 2025:
Tasklist: 
- Create a script to BLAST ORIs --> done, 134 ORIs identified using BLAST
- Append AMR info to df
- Figures to generate:
  - Length distribution
  - Circle plot of association for ORIs


# Jul 19, 2025:
- ORI stats:
  - 136 plasmids have a type 1 ORI, 431 with type 1 or 2 ORI
- Create a script to find the cooccurrence frequency of identified ORIs and how often each one occurs


# Jul 21, 2025:
- Generating a length distribution plot for 1000 plasmids, multiple ORIs vs TOTOs, includes type 1 and 2 ORIs
- empty ORIs: 997, empty BLASTs: 33374
- Running AMRfinder on entire DB (75000), running heterodimer-finder.py on entire DB
  - Need to figure out a way to append this data into a final df that has all the relevant info necessary
- Launch a jupyter session from vscode: `sbatch vscode.sh` `tail -f slurm-7885318.out`
- Q: How do we know these are heterodimers and not homodimers (same ORI sandwiched between IS)


# Jul 22, 2025:
- Figure out why heterodimer-finder.py isn't working for full DB
  - Seems final ORI files are empty after parsing orivfinder outputs --> rerunning origins.py
    - Problem was that some of the orivfinder outputs are tsvs and some are csvs, also for the tsvs some columns had commas which was messing up the line reading... a whole mess so some files have been skipped
- Generate plots
  - Generated a length distribution plot, realized sub dfs are lackluster because original tsv file has multiple ids for the same plasmid which is screwing up the format of the data --> fixed this, new tsv called `final_df.tsv`, append data to this
- Append data and make a final df with all info
- Reran for empty files: empty ORIs: 12383, empty BLASTs: 33374

To Do:
- Append AMR info to final dataframe
- I think I deleted the final/ORI files on accident (judging from the ORI empty df concat error when running heterodimer-finder.py)-- run origins.py then how_many.py then heterodimer-finder.py again...
  - Also possibly next time sbatch heterodimer-finder.py because it takes a while to run lol
  - NVM IT RAN!!!! plasmid hits: 9817, total hits: 21674
  - mORIs: 18350 distinct plasmids (type 1 ORIs only)


# Jul 23, 2025:
- ORI stats:
  - 37963 plasmids with a type 1 ORI, of those-- 18350 with 2+ ORIs and 9817 TOTOs
- Length stats:
  - avg all: 92127.55898123325 --> bimodal distribution, find the peak of each peak
  - avg mORIs: 132316.43543543544
  - avg TOTOs: 177994.69
- Some statistical analyses:
  - for lengths-- randomly sample pairs of 10,000 for which one is longer. This will generate a p-value
    - pretty high p-val --> not surprising as this is suggestive of TOTOs being a subset of mORIs which is true. also no significant difference (in fact is higher for exclusive which is contrary to intuition) between mORIs inclusive and exclusive of TOTOs
  - another metric is to bin by the number of ORIs --> this makes sense because intuitively TOTOs should have more ORIs since they can be the fusion of plasmids that originally had several ORIs 
    - can do the same for # of ins seqs on a plasmid
  - for AMRs-- randomly sample 10,000 plasmids with replacement, repeat 100 trials, keep track of the fraction with AMRs and plot between different groups
  - Want to be able to prove that the output results (TOTOs) are a real result and not something due to biases --> so remove all the possible biases


# Jul 24, 2025:
- Only have 134 unique type 1 ORIs after using BLAST to identify --> extend to type 2 ORIs as well
- Cooccurrence table (including type 2 ORIs)
  - Possibly rerun heterodimer-finder.py inclusive of type 2 ORIs
- TOOTs (as a direct control) and type 2 ORIs
- Work on slides


# Jul 25, 2025:
- BLAST ORIs once that finishes running
- Interesting Qs:
  - Which IS and ORI are seen together
  - Compare the expected and actual for co-occurence freqs


# Jul 26, 2025:
- BLASTed ORIs for type 1,2,3. Parsed output and assigned highest hits
- Created script to find the absolute freq of each ORI, script for finding co-occurence matrix in progress
  - When calling this func, append type 1 and type 2 outputs together
- Ran heterodimer-finder.py on type 2 ORIs, analyze this data and maybe put together with type 1

# Jul 28, 2025:
- Split cooccurence into batch jobs, once they finish running stitch all the matrices together
- Most important: work on slides


# Jul 29, 2025:
- BLAST ORI identification script is wrong... fix this bc each plasmid is only associated with 1 ORI rn
  - Don't know where the bug is :( go through each step or maybe start from scratch


# Aug 4, 2025:
- For identifying ORIs, use the mmseqs clustering approach (maybe back to 90% identity, 80% coverage) and BLAST the representative seqs of the clusters against OriVFinder DB
  - BLASTing ORI seqs against OriVFinder DB did not work as there are many more de novo predictions that expected.


# Aug 28, 2025:
- Can't access O2 portal --> access revoked?
- Figure out where I left on the cooccurence table and grouping/clustering origins


# Aug 29, 2025:
- Cluster ORIs with mmseqs2 and see which clusters are identified using BLAST
  - if this works do a quick analysis of co-occurence using pipeline
- Email Marisela and ask abt credit
