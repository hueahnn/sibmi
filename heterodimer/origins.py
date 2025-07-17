# author: hueahnn
# begin: 07/02/2025
# updated: 07/17/25
# purpose: determine how many plasmids have multiple origins and the co-occurence of origins

import pandas as pd
import os
import sys

def main(PLASMID_PATH):
    OUTPUT_PATH = f"ORIs.summary.tsv"
    COMBINED_PATH = f"ORIs.df.tsv"
    PLASMIDS = []
    with open(PLASMID_PATH, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]

    # determine how many plasmids have multiple origins
    counter = 0 # counter for how many plasmids have more than one origin
    combined_df = pd.DataFrame(columns=["seqID", "Type", "Intergenic_Start", "Intergenic_End", "Intergenic_Sequence"])
    for PLASMID in PLASMIDS:
        ORI_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
        if os.path.getsize(ORI_PATH) != 0:
            df = pd.read_csv(ORI_PATH, sep="\t")
            df = df[["seqID", "Type", "Intergenic_Start", "Intergenic_End", "Intergenic_Sequence"]]
            if (len(df) > 1):
                counter+=1
                combined_df = pd.concat([df, combined_df])
                with open(OUTPUT_PATH, "a") as f:
                    print(PLASMID,file=f)
                    print(df, file=f)
    combined_df.to_csv(COMBINED_PATH, sep="\t", index=False)
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

# convert ORI df to fasta file for clustering with mmseqs2
def mmseqs(INPUT_FILE, OUTPUT_FILE):
    df = pd.read_csv(INPUT_FILE, sep="\t")
    count = 0
    with open(OUTPUT_FILE, "w") as f:
        for _, row in df.iterrows():
            header = f">{row['seqID']}_{count}"
            seq = row["Intergenic_Sequence"]
            f.write(f"{header}\n{seq}\n")
            count+=1

            
if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])





