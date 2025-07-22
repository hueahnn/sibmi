# author: hueahnn
# begin: 07/11/25
# updated 07/11/25
# purpose: how many files are empty?

import pandas as pd
import sys
import os



def main(TXT_FILE):
    PLASMIDS = []
    with open(TXT_FILE, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    df = pd.read_csv("../filtered.ALL.plasmid_list.download.tsv", sep="\t")
    df = df[df.Plasmid_ID.isin(PLASMIDS)]
    df.to_csv("test.csv", sep="\t")
    missing = ~df.Plasmid_ID.isin(PLASMIDS)
    missing.to_csv("missing.csv", sep="\t")
    BLAST_FILE = "empty_BLASTs.txt"
    open(BLAST_FILE, "w").close()
    empty_BLASTs = 0
    ORI_FILE = "empty_ORIs.txt"
    open(ORI_FILE, "w").close()
    empty_ORIs = 0
    for plasmid in PLASMIDS:
        ORI_PATH = f"heterodimer/ORIs/{plasmid}/All_IGSs.csv"
        open(ORI_PATH, "a").close()
        BLAST_PATH = f"heterodimer/{plasmid}.pairwise.blast.tsv"
        open(BLAST_PATH, "a").close()
        if (os.path.getsize(BLAST_PATH) == 0):
            empty_BLASTs += 1
            with open(BLAST_FILE, "a") as f:
                print(plasmid, file=f)
        if (os.path.getsize(ORI_PATH) == 0):
            empty_ORIs += 1
            with open(ORI_FILE, "a") as f:
                print(plasmid, file=f)
        else:
            ori_df = pd.read_csv(ORI_PATH, sep="\t")
            if ori_df.empty:
                empty_ORIs += 1
                with open(ORI_FILE, "a") as f:
                    print(plasmid, file=f)
            elif ori_df.iloc[0,0] == plasmid:
                empty_ORIs += 1
                with open(ORI_FILE, "a") as f:
                    print(plasmid, file=f)
                open(ORI_PATH, "w").close()
    print(f"empty ORIs: {empty_ORIs}\nempty BLASTs: {empty_BLASTs}")



if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
