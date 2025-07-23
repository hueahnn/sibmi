# author: hueahnn
# begin: 07/02/2025
# updated: 07/23/25
# purpose: determine how many plasmids have multiple origins and the co-occurence of origins

import pandas as pd
import os
import sys
import csv

def multiple(PLASMID_PATH):
    OUTPUT_PATH = f"mORIs.summary.tsv"
    COMBINED_PATH = f"mORIs.df.tsv"
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


# reading in a mix of tsvs and csvs
def read_auto_delim_csv(filepath):
    import csv
    with open(filepath, 'r', newline='') as f:
        sample = f.read(2048)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t', ';'])
            delimiter = dialect.delimiter
        except csv.Error:
            delimiter = ','  # fallback
    # use python engine and warn if rows are skipped
    try:
        df = pd.read_csv(filepath, delimiter=delimiter, engine='python', on_bad_lines='warn')
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
        return pd.DataFrame()
    return df

def cleanup_orivfinder(PLASMID):
    OUTPUT_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
    open(OUTPUT_PATH, "w").close()
    PATH = f"heterodimer/ORIs/{PLASMID}/"
    input_file = os.path.join(PATH, "All_IGSs.csv")
    if (os.path.getsize(f"{PATH}/All_IGSs.csv") == 0):
        return
    igs_df = read_auto_delim_csv(input_file)
    if igs_df.empty or "Type" not in igs_df.columns:
        return
    # clean up and stitch together igs_df
    igs_df["seqID"] = PLASMID
    igs_df = igs_df[igs_df.Type==1]
    if "Unnamed: 0" in igs_df.columns:
        igs_df = igs_df.drop(columns=["Unnamed: 0"])
    igs_df.to_csv(OUTPUT_PATH, sep="\t", index=False)


def cleanup_orivfinder_multiple(FILE):
    PLASMIDS = []
    with open(FILE, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    for PLASMID in PLASMIDS:
        OUTPUT_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
        open(OUTPUT_PATH, "w").close()
        PATH = f"heterodimer/ORIs/{PLASMID}/"
        input_file = os.path.join(PATH, "All_IGSs.csv")
        if (os.path.getsize(f"{PATH}/All_IGSs.csv") == 0):
            continue
        igs_df = read_auto_delim_csv(input_file)
        if igs_df.empty or "Type" not in igs_df.columns:
            continue
        # clean up and stitch together igs_df
        igs_df["seqID"] = PLASMID
        igs_df = igs_df[igs_df.Type==1]
        if "Unnamed: 0" in igs_df.columns:
            igs_df = igs_df.drop(columns=["Unnamed: 0"])
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


# assign highest BLAST hit
def identity(FILE):
    OUTPUT_FILE = "type.one.ORIs.identified.tsv"
    df = pd.read_csv(FILE, sep="\t")
    df = df.loc[df.groupby("qseqid")["bitscore"].idxmax()]
    df.to_csv(OUTPUT_FILE, sep="\t")


# create empty file if missing
def missing(FILE):
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    for plasmid in PLASMIDS:
        ORI_PATH = f"heterodimer/final/{plasmid}.All_IGSs.csv"
        dir_path = os.path.dirname(ORI_PATH)
        os.makedirs(dir_path, exist_ok=True)
        open(ORI_PATH, 'a').close()

def num_origins(FILE):
    ORI_df = pd.read_csv("ORI.summary.tsv", sep="\t", index_col=False)
    main_df = pd.read_csv(FILE, sep="\t")
    main_df["num_ORIs"] = 0
    for id in ORI_df["seqID"].unique():
        sub_df = ORI_df[ORI_df["seqID"]==id]
        ORIs = len(sub_df)
        main_df.loc[main_df["Plasmid_ID"]==id, "num_ORIs"] = ORIs
    main_df.to_csv(FILE, sep="\t", index=False)




if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])





