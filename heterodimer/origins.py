# author: hueahnn
# begin: 07/02/2025
# updated: 08/04/25
# purpose: determine how many plasmids have multiple origins and the co-occurence of origins

import pandas as pd
import os
import sys
import csv
from tqdm import tqdm
from Bio import SeqIO
from collections import Counter
import numpy as np
from itertools import combinations_with_replacement



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
        OUTPUT_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.3.csv"
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
        igs_df = igs_df[igs_df.Type==3]
        if "Unnamed: 0" in igs_df.columns:
            igs_df = igs_df.drop(columns=["Unnamed: 0"])
        igs_df.to_csv(OUTPUT_PATH, sep="\t", index=False)


# convert ORI df to fasta file
def fasta(INPUT_FILE, OUTPUT_FILE):
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

    df = pd.read_csv(FILE, sep="\t")
    for id in tqdm(range(len(df["qseqid"].unique()))):
        sub = df[df["qseqid"]==id]
        sub = sub.sort_values("bitscore", ascending=False)
        
    df.to_csv(FILE, sep="\t")


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

# how many ORIs does each plasmid have?
def num_origins(FILE):
    ORI_df = pd.read_csv("ORI.summary.tsv", sep="\t", index_col=False)
    main_df = pd.read_csv(FILE, sep="\t")
    main_df["num_ORIs"] = 0
    for id in ORI_df["seqID"].unique():
        sub_df = ORI_df[ORI_df["seqID"]==id]
        ORIs = len(sub_df)
        main_df.loc[main_df["Plasmid_ID"]==id, "num_ORIs"] = ORIs
    main_df.to_csv(FILE, sep="\t", index=False)

def num_origins(FILE, DF):
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    main_df = pd.read_csv(DF, sep="\t")
    main_df["num_ORIs"] = 0
    for plasmid in PLASMIDS:
        PATH = f"heterodimer/final/{plasmid}.All.IGSs.csv"
        if os.path.getsize(PATH) == 0:
            continue
        sub_df = pd.read_csv(PATH, sep="\t")
        ORIs = len(sub_df)
        main_df.loc[main_df["Plasmid_ID"]==plasmid, "num_ORIs"] = ORIs
    main_df.to_csv(DF, sep="\t", index=False)


# join all orivfinder outputs into one df
def ori_df(FILE):
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    OUTPUT_FILE = "all.ORIs.tsv"
    df = pd.DataFrame(columns=[])
    for plasmid in tqdm(PLASMIDS):
        PATH = f"heterodimer/ORIs/{plasmid}/All_IGSs.csv"
        if os.path.getsize(PATH) == 0:
            continue
        plasmid_df = pd.read_csv(PATH, sep="\t",index_col=False)
        df = pd.concat([df, plasmid_df], ignore_index=True)
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)


# joining all the final ORI data into one df
def origin_df(FILE):
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    OUTPUT_FILE = "ORIs.tsv"
    ori_list = []
    for plasmid in tqdm(PLASMIDS):
        PATH1 = f"heterodimer/final/{plasmid}.All.IGSs.csv"
        PATH2 = f"heterodimer/final/{plasmid}.All.IGSs.2.csv"
        PATH3 = f"heterodimer/final/{plasmid}.All.IGSs.3.csv"
        if os.path.getsize(PATH1) == 0:
            continue
        if os.path.getsize(PATH2) == 0:
            continue
        if os.path.getsize(PATH3) == 0:
            continue
        sub1_df = pd.read_csv(PATH1, sep="\t")
        sub2_df = pd.read_csv(PATH2, sep="\t")
        sub3_df = pd.read_csv(PATH3, sep="\t")
        ori_list.append(sub1_df)
        ori_list.append(sub2_df)
        ori_list.append(sub3_df)
    ori_df = pd.concat(ori_list, ignore_index=True)
    ori_df.to_csv(OUTPUT_FILE, sep="\t")
        

# finding the abs freq of each ORI
def abs_freq(FILE):
    # list of every ORI in OriVFinder DB
    ORIs = list(set([record.id for record in SeqIO.parse("../orivfinder/app/OriVDB/OriV.fasta", "fasta")]))
    # find the absolute freq of each ORI
    ABS_OUTPUT_FILE = "ORIs.abs.freq.tsv"
    abs_df = pd.DataFrame(ORIs, columns=["ORI"])
    abs_df["count"] = 0
    identified_df = pd.read_csv(FILE, sep="\t")
    identified_df = identified_df.drop_duplicates()
    identified_df["id"] = identified_df["qseqid"].str.split(".").str[1]
    identified_ORIs = identified_df["sseqid"]
    counts = Counter(identified_ORIs)
    for ori in tqdm(ORIs):
        abs_df.loc[abs_df["ORI"]==ori, "count"] = counts[ori]
    abs_df.to_csv(ABS_OUTPUT_FILE, sep="\t")


# building a co-occurrence matrix for all ORIs, arg: pass in csv of identified ORIs
def cooccurrence(FILE, COOCCUR_OUTPUT_FILE):
    # list of every ORI in OriVFinder DB
    ORIs = list(set([record.id for record in SeqIO.parse("../orivfinder/app/OriVDB/OriV.fasta", "fasta")]))
    identified_df = pd.read_csv("type.1and2.ORIs.identified.tsv", sep="\t")
    identified_df["id"] = identified_df["qseqid"].str.split(".").str[0]
    # co-occurence table
    matrix = np.zeros((len(ORIs), len(ORIs)))
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    for i in tqdm(range(len(ORIs))):
        for j in tqdm(range(i,len(ORIs))):
            count = 0
            for plasmid in tqdm(PLASMIDS):
                sub_df = identified_df[identified_df["id"]==plasmid]
                plasmid_oris = sub_df["sseqid"]
                if (ORIs[i] in plasmid_oris) and (ORIs[j] in plasmid_oris):
                    count+=1
                    continue
            matrix[i,j] = count
    matrix_df = pd.DataFrame(matrix)
    matrix_df.to_csv(COOCCUR_OUTPUT_FILE, sep="\t")


# a better implementation for finding the cooccurrences
def cooccurrence2(FILE, OUTPUT_FILE):
    # list of every ORI in OriVFinder DB
    ORIs = list(set([record.id for record in SeqIO.parse("../orivfinder/app/OriVDB/OriV.fasta", "fasta")]))
    ORI_dict = {ori: i for i, ori in enumerate(ORIs)}
    identified_df = pd.read_csv("type.1and2.ORIs.identified.tsv", sep="\t")
    identified_df["id"] = identified_df["qseqid"].str.split("_").str[:2].str.join("_")
    print(identified_df)
    # co-occurence table
    matrix = np.zeros((len(ORIs), len(ORIs)))
    PLASMIDS = []
    with open(FILE, 'r') as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    print(PLASMIDS)
    for plasmid in tqdm(PLASMIDS):
        sub_df = identified_df[identified_df["id"]==plasmid]
        print(sub_df)
        plasmid_oris = sub_df["sseqid"]
        pairs = list(combinations_with_replacement(plasmid_oris, 2))
        for pair in tqdm(pairs):
            print(pair)
            ori1 = ORI_dict[pair[0]]
            ori2 = ORI_dict[pair[1]]
            matrix[ori1,ori2]+=1
    matrix_df = pd.DataFrame(matrix)
    matrix_df.to_csv(OUTPUT_FILE, sep="\t")






if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])