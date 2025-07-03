# author: hueahnn
# begin: 07/02/2025
# purpose: determine how many plasmids have multiple origins and the co-occurence of origins

import pandas as pd
import os
import sys

def main(PLASMID_PATH, FILE_PATH):
    PLASMIDS = []
    OUTPUT_PATH = f"{FILE_PATH}/ORIs.summary.tsv"
    with open(PLASMID_PATH, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]

    # determine how many plasmids have multiple origins
    counter = 0 # counter for how many plasmids have more than one origin
    # combined_df = pd.DataFrame("SEQUENCE", "PRODUCT")
    for plasmid in PLASMIDS:
        df = pd.read_csv(f"{FILE_PATH}/{plasmid}.ORIs.tsv", sep="\t")
        df = df[["SEQUENCE", "PRODUCT"]]
        if (len(df) > 1):
            counter+=1
            with open(OUTPUT_PATH, "a") as f:
                print(plasmid,file=f)
                print(df, file=f)
            # combined_df = pd.concat([df, combined_df])
    # combined_df.to_csv(OUTPUT_PATH, sep="\t", index=False)
    with open(OUTPUT_PATH, "a") as f:
        print(counter, file=f)

            
if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])





