# author: hueahnn
# begin: 06/24/25
# purpose: for filtering for the heterodimers of interest 
import pandas as pd
import os
import sys

### implementation for running on a single plasmid (for snakemake) ###
def single(path, plasmid):
    # import .tsv files containing ORIs and IS to dataframes]
    ins_seq_filepath = f"{path}/{plasmid}/fasta/{plasmid}.fasta.tsv"
    ori_filepath = f"{path}/ORIs/{plasmid}.ORIs.tsv"
    open(f"{path}/{plasmid}-heterodimer-outputs.txt", "w").close()
    if (os.path.getsize(ins_seq_filepath) != 0):
        ins_seq_df = pd.read_csv(ins_seq_filepath, sep='\t')
        ins_seq_df = ins_seq_df.rename(columns={"isBegin":"begin", "isEnd":"end", "cluster":"feature"})
        ins_seq_df["type"] = "is"
        ori_df = pd.read_csv(ori_filepath, sep='\t')
        ori_df = ori_df.rename(columns={"SEQUENCE":"seqID", "START":"begin", "END":"end", "GENE":"feature"})
        ori_df["type"] = "ori"
        # specify output file directory
        output_file = f"{path}/{plasmid}-heterodimer-outputs.txt"
        count(ins_seq_df, ori_df, output_file)    
        if (os.path.getsize(output_file) != 0):
            with open(f"{path}/heterodimer-outputs.txt", "a") as f:
                print(plasmid, file=f)

# counting the number of class changes for each insertion sequence within a plasmid
def count(ins_seq_df, ori_df, output_file):
    # filter out seqs with 1 IS and combine info into one df
    ins_seq_df = ins_seq_df[ins_seq_df.ncopy4is != 1]

    # rename cols to match and combine the dfs into one df
    ori_simple = ori_df[["seqID", "begin", "end", "type", "feature"]]
    is_simple = ins_seq_df[["seqID", "begin", "end", "type", "feature"]]
    df = pd.concat([is_simple, ori_simple])

    # determine sequential order of features + number of class changes + output seqID
    for id in df["seqID"].unique():
        plasmid_df = df[df.seqID == id]
        plasmid_df = plasmid_df.sort_values("begin")
        is_only = plasmid_df[plasmid_df.type == "is"] 
        for feat in is_only["feature"].unique():
            is_df = plasmid_df[(plasmid_df.type == "ori") | ((plasmid_df.feature == feat) & (plasmid_df.type == "is"))]
            changes = 0
            for i in range(0,len(is_df)):
                curr = is_df.type.iat[i]
                next = None
                if (i == len(is_df)-1):
                    next = is_df.type.iat[0]
                else:
                    next = is_df.type.iat[i+1]
                if (curr != next) :
                    changes += 1
            if (changes > 2):
                with open(output_file, "a") as f:
                    print(f"{id}\t{feat}\t{changes}", file=f)
                    print(is_df, file=f)


### below is an implementation for running the script on one giant aggregated file containing the info for multiple plasmids instead of a single plasmid ###
def aggregated(): 
    # import .tsv files containing ORIs and IS to dataframes
    path = [f for f in os.listdir("test-batch-1/ins-seqs/fusion-seqs") if f.endswith(".tsv")]
    ins_seq_df = pd.concat((pd.read_csv(f"test-batch-1/ins-seqs/fusion-seqs/{file}", sep='\t') for file in path), ignore_index=True)
    ins_seq_df = ins_seq_df.rename(columns={"isBegin":"begin", "isEnd":"end", "cluster":"feature"})
    ins_seq_df["type"] = "is"
    # todo: later on add a line to export the df into a tsv
    ori_df = pd.read_csv("test-batch-1/ORIs.tsv", sep='\t')
    ori_df = ori_df.rename(columns={"SEQUENCE":"seqID", "START":"begin", "END":"end", "GENE":"feature"})
    ori_df["type"] = "ori"
    # specify output file directory
    output_file = "test-batch-1/heterodimer-outputs.txt"
    count(ins_seq_df, ori_df, output_file)

    

if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
