# author: hueahnn
# begin: 06/30/2025
# purpose: parse BLAST output table + run a pairwise BLAST 
import pandas as pd
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import string

def main(PLASMID):
    # for each IS found, determine if they are in overlapping regions and for ones that are not run a pairwise BLAST to determine if they are the same
        OUTPUT_FILE = f"heterodimer/{PLASMID}.unique.hits.tsv"
        open(OUTPUT_FILE, "w").close()
        PAIRWISE_OUTPUT_FILE = f"heterodimer/{PLASMID}.pairwise.blast.tsv"
        open(PAIRWISE_OUTPUT_FILE, "w").close()
        # for files with no hits
        BLAST_FILEPATH = f"heterodimer/{PLASMID}.blast.tsv"
        if (os.path.getsize(BLAST_FILEPATH) == 0):
            return
        df = pd.read_csv(f"heterodimer/{PLASMID}.blast.tsv", sep="\t", header=None)
        df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        df = df.sort_values("bitscore", ascending=False)
        # df.to_csv("test.blast.tsv", sep="\t")
        unique = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

        # step 1: some regions have several hits-- filter these out by checking the position of the region and only keep the highest bitscore hit
        skip = []
        for i in range(0,len(df)):
            if i in skip:
                continue
            else:
                temp = df.iloc[i]
                unique.loc[len(unique)] = temp
                start = df.qstart.iat[i]
                end = df.qend.iat[i]
                for j in range(i,len(df)):
                    if (df.qend.iat[j] > start) & (df.qstart.iat[j] < end):
                        skip.append(j)
        unique.to_csv(OUTPUT_FILE, sep="\t")
        
        # step 2: for every unique insertion sequence, run a pairwise BLAST with all ins seq to determine if they actually are the same ins seq (cutoff 95%)
        db = PLASMID.split('_')[0]
        fasta = f"../fasta-files/{db}/fasta/{PLASMID}.fasta"
        record = next(SeqIO.parse(fasta, "fasta"))
        sequence = str(record.seq)
        # generate fasta file of hits
        HITS_FILE = f"heterodimer/{PLASMID}.hits.fasta"
        records = []
        for i in range(0,len(unique)):
            start = unique.qstart.iat[i]
            end = unique.qend.iat[i]
            ins_seq = sequence[min(start,end): max(start,end)]
            record = SeqRecord(Seq(ins_seq), id=f"{unique.sseqid.iat[i]}_{i}", description="")
            records.append(record)
        SeqIO.write(records, HITS_FILE, "fasta")
        # upload as db and run BLAST
        os.system(f"makeblastdb -in {HITS_FILE} -dbtype nucl -out blast-dbs/{PLASMID}")
        os.system(f"blastn -db blast-dbs/{PLASMID} -query {HITS_FILE} -outfmt 6 >> {PAIRWISE_OUTPUT_FILE}")
        # reformat output
        if (os.path.getsize(PAIRWISE_OUTPUT_FILE) != 0):
            reformat = pd.read_csv(PAIRWISE_OUTPUT_FILE, sep="\t", header=None)
            reformat.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            reformat.to_csv(PAIRWISE_OUTPUT_FILE, sep="\t")

def parse(PLASMID):
    # for all the pident 100, add unique ids to list
    PAIRWISE_INPUT_FILE = f"heterodimer/{PLASMID}.pairwise.blast.tsv"
    UNIQUE_INPUT_FILE = f"heterodimer/{PLASMID}.unique.hits.tsv"
    OUTPUT_FILE = f"heterodimer/final/{PLASMID}.blast.tsv"
    open(OUTPUT_FILE, "w").close()
    if (os.path.getsize(PAIRWISE_INPUT_FILE) == 0):
            return
    pairwise_df = pd.read_csv(PAIRWISE_INPUT_FILE, sep="\t")
    unique_df = pd.read_csv(UNIQUE_INPUT_FILE, sep="\t")
    if "Unnamed: 0" in unique_df.columns:
        unique_df = unique_df.drop(columns=["Unnamed: 0"])
    groups = []
    g = set()
    for i in range(0,len(unique_df)):
         same = []
         s = set()
         for j in range(0,len(pairwise_df)):
            qid = pairwise_df.qseqid.iat[j].split("_")
            qid = int(qid[len(qid)-1])
            sid = pairwise_df.sseqid.iat[j].split("_")
            sid = int(sid[len(sid)-1])
            if i not in s:
                same.append(i)
                s.add(i)
            if (qid==i) & (sid!=i):
                same.append(sid)
         same = sorted(same)
         if tuple(same) not in g:
            groups.append(same)
            g.add(tuple(same))
    append_group = [None]*len(unique_df)
    gen = 0
    for i in range(0,len(groups)):
         for num in groups[i]:
              append_group[num] = f"{PLASMID}_{gen}"
         gen += 1
    unique_df["group"] = append_group
    unique_df.to_csv(OUTPUT_FILE,sep="\t",index=False)

def parse_multiple(FILE):
    PLASMIDS = []
    with open(FILE, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    for PLASMID in PLASMIDS:
        PAIRWISE_INPUT_FILE = f"heterodimer/{PLASMID}.pairwise.blast.tsv"
        UNIQUE_INPUT_FILE = f"heterodimer/{PLASMID}.unique.hits.tsv"
        OUTPUT_FILE = f"heterodimer/final/{PLASMID}.blast.tsv"
        open(OUTPUT_FILE, "w").close()
        if (os.path.getsize(PAIRWISE_INPUT_FILE) == 0):
            continue
        pairwise_df = pd.read_csv(PAIRWISE_INPUT_FILE, sep="\t")
        unique_df = pd.read_csv(UNIQUE_INPUT_FILE, sep="\t")
        if "Unnamed: 0" in unique_df.columns:
            unique_df = unique_df.drop(columns=["Unnamed: 0"])
        groups = []
        g = set()
        for i in range(0,len(unique_df)):
            same = []
            s = set()
            for j in range(0,len(pairwise_df)):
                qid = pairwise_df.qseqid.iat[j].split("_")
                qid = int(qid[len(qid)-1])
                sid = pairwise_df.sseqid.iat[j].split("_")
                sid = int(sid[len(sid)-1])
                if i not in s:
                    same.append(i)
                    s.add(i)
                if (qid==i) & (sid!=i):
                    same.append(sid)
            same = sorted(same)
            if tuple(same) not in g:
                groups.append(same)
                g.add(tuple(same))
        append_group = [None]*len(unique_df)
        gen = 0
        for i in range(0,len(groups)):
         for num in groups[i]:
              append_group[num] = f"{PLASMID}_{gen}"
         gen += 1
        unique_df["group"] = append_group
        unique_df.to_csv(OUTPUT_FILE,sep="\t",index=False)



if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
