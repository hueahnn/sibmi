# author: hueahnn
# begin: 07/02/2025
# updated: 07/03/25
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

def cleanup_orivfinder(PLASMID):
    OUTPUT_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
    open(OUTPUT_PATH, "w").close()
    PATH = f"heterodimer/ORIs/{PLASMID}/"
    if (os.path.getsize(f"{PATH}/All_IGSs.csv") == 0):
        return
    igs_df = pd.read_csv(f"{PATH}/All_IGSs.csv", sep="\t")
    igs_df["seqID"] = PLASMID
    rip_df = pd.read_csv(f"{PATH}/RIP.csv", sep="\t")
    # clean up and stitch together igs_df
    if igs_df.empty:
        return
    igs_df = igs_df[igs_df.Type < 3]
    if "Unnamed: 0" in igs_df.columns:
        igs_df = igs_df.drop(columns=["Unnamed: 0"])
    igs_df.to_csv(OUTPUT_PATH, sep="\t",index=False)


def cleanup_orivfinder_multiple(FILE):
    PLASMIDS = []
    with open(FILE, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    for PLASMID in PLASMIDS:
        OUTPUT_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
        open(OUTPUT_PATH, "w").close()
        PATH = f"heterodimer/ORIs/{PLASMID}/"
        if (os.path.getsize(f"{PATH}/All_IGSs.csv") == 0):
            continue
        igs_df = pd.read_csv(f"{PATH}/All_IGSs.csv", sep="\t")
        igs_df["seqID"] = PLASMID
        rip_df = pd.read_csv(f"{PATH}/RIP.csv", sep="\t")
        # clean up and stitch together igs_df
        if igs_df.empty:
            continue
        igs_df = igs_df[igs_df.Type < 3]
        if "Unnamed: 0" in igs_df.columns:
            igs_df = igs_df.drop(columns=["Unnamed: 0"])
        # columns to append
        igs_df.to_csv(OUTPUT_PATH, sep="\t", index=False)



            
if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])





