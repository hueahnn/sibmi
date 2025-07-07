# author: hueahnn
# begin: 06/30/2025
# purpose: parse BLAST output table + run a pairwise BLAST 
import pandas as pd
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main(FILE_PATH):
    #PLASMIDS = []
    #with open(PLASMID_PATH, "r") as f:
    #    PLASMIDS = [line.strip() for line in f if line.strip()]

    # for each IS found, determine if they are in overlapping regions and for ones that are not run a pairwise BLAST to determine if they are the same
    PLASMIDS = ["DDBJ_AP017321.1"]
    for plasmid in PLASMIDS:
        df = pd.read_csv(f"{FILE_PATH}/{plasmid}.blast.tsv", sep="\t", header=None)
        df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        df = df[df.bitscore > 100]
        df = df.sort_values("bitscore", ascending=False)
        # df.to_csv("test.blast.tsv", sep="\t")
        unique = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

        # step 1: some regions have several hits -- filter out these by checking the position of the region and only keep the highest bitscore hit
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
        unique.to_csv("test.blast.tsv", sep="\t")
        #with open("test.blast.tsv", "a") as f:
            #print(unique, file=f)
        

        # step 2: for every unique insertion sequence, run a pairwise BLAST on each pair to determine if they actually are the same ins seq (cutoff 90% or so)
        plasmid = unique.qseqid.iat[0]
        db = plasmid.split('_')[0]
        PAIRWISE_OUTPUT_FILE = f"heterodimer/{plasmid}.pairwise.blast.tsv"
        fasta = f"../fasta-files/{db}/fasta/{plasmid}.fasta"
        open(PAIRWISE_OUTPUT_FILE, "w").close()
        # upload plasmid seq as db
        os.system(f"makeblastdb -in {fasta} -dbtype nucl -out blast-dbs/{plasmid}")
        record = next(SeqIO.parse(fasta, "fasta"))
        sequence = str(record.seq)
        for i in range(0,len(unique)):
            start = unique.qstart.iat[i]
            end = unique.qend.iat[i]
            ins_seq = sequence[min(start,end): max(start,end)]
            # write to temp query file
            temp_query_file = "temp_query.fasta"
            SeqIO.write(SeqRecord(Seq(ins_seq), id=f"{plasmid}_{i}", description=""), temp_query_file, "fasta")
            # run BLAST
            os.system(f"blastn -db blast-dbs/{plasmid} -query {temp_query_file} -outfmt 6 >> {PAIRWISE_OUTPUT_FILE}")
        reformat = pd.read_csv(PAIRWISE_OUTPUT_FILE, sep="\t", header=None)
        reformat.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        reformat = reformat[reformat.pident >= 90]
        reformat.to_csv(PAIRWISE_OUTPUT_FILE, sep="\t")




if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
