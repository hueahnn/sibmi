# author: hueahnn
# begin: 06/30/2025
# purpose: parse BLAST output table + run a pairwise BLAST 
import pandas as pd
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
        reformat = pd.read_csv(PAIRWISE_OUTPUT_FILE, sep="\t", header=None)
        reformat.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        reformat.to_csv(PAIRWISE_OUTPUT_FILE, sep="\t")


if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
