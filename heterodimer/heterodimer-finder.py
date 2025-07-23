# author: hueahnn
# begin: 06/24/25
# updated: 07/23/25
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
def count(ins_seq_df, ori_df, output):
    output_file = f"heterodimer-outputs-{output}.txt"
    open(output_file, "w").close()
    output_df = pd.DataFrame(columns=["seqID", "begin", "end", "type", "feature"])
    output_plasmids = f"heterodimer-plasmids-{output}.txt"
    # filter out seqs with 1 IS and combine info into one df
    ins_seq_df = ins_seq_df[ins_seq_df.copies > 1]

    # rename cols to match and combine the dfs into one df
    ori_simple = ori_df[["seqID", "begin", "end", "type", "feature"]]
    is_simple = ins_seq_df[["seqID", "begin", "end", "type", "feature"]]
    df = pd.concat([is_simple, ori_simple])

    # determine sequential order of features + number of class changes + output seqID
    total = 0
    plasmids = 0
    for id in df["seqID"].unique():
        plasmid_df = df[df.seqID == id]
        plasmid_df = plasmid_df.sort_values("begin")
        is_only = plasmid_df[plasmid_df.type == "is"] 
        first = 0
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
                total+=1
                if first==0:
                    plasmids+=1
                    with open(output_plasmids, "a") as f:
                        print(id, file=f)
                first+=1
                with open(output_file, "a") as f:
                    print(f"{id}\t{feat}\t{changes}", file=f)
                    print(is_df, file=f)
                output_df = pd.concat([output_df,is_df], ignore_index=True)
    with open(output_file, "a") as f:
        print(f"plasmid hits: {plasmids}\ntotal hits: {total}",file=f)
    output_df.to_csv(f"TOTO_df_{output}.csv",sep="\t")


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


### for input txt file containing ids of every seq
def aggregated(FILE, OUTPUT_NUM):
    PLASMIDS = []
    with open(FILE, "r") as f:
        PLASMIDS = [line.strip() for line in f if line.strip()]
    ins_seq_df = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'group', 'copies'])
    ori_df = pd.DataFrame(columns=['Accession_Number', 'Evidence', 'Type', 'Intergenic_Start', 'Intergenic_End', 'Intergenic_Sequence', 'Iteron_Entropy_Score', 'Iteron_df', 'Pattern_Score', 'Pattern_Df', 'AT_Score', 'AT_Df', 'RNA_df', 'Total_Score', 'Sum_Score', 'gene', 'gene_id', 'gene_start', 'gene_end', 'mmseqs_hit', 'seqID'])
    for PLASMID in PLASMIDS:
        IS_PATH = f"heterodimer/final/{PLASMID}.blast.tsv"
        ORI_PATH = f"heterodimer/final/{PLASMID}.All.IGSs.csv"
        if (os.path.getsize(IS_PATH) != 0) and (os.path.getsize(ORI_PATH) != 0):
            ins_seq_df = pd.concat([ins_seq_df,pd.read_csv(IS_PATH, sep='\t')], ignore_index=True)
            ori_df = pd.concat([ori_df,pd.read_csv(ORI_PATH, sep='\t')], ignore_index=True)
    # put together ins seq data
    ins_seq_df = ins_seq_df.rename(columns={"qseqid":"seqID","qstart":"begin", "qend":"end", "group":"feature"})
    ins_seq_df["type"] = "is"
    # do some filtering for eval, length, percent identity, and such?
    IS_OUTPUT = "IS.summary.tsv"
    ins_seq_df.to_csv(IS_OUTPUT, sep="\t", index=False)
    # put together ORI data
    ori_df = ori_df.rename(columns={"Intergenic_Start":"begin", "Intergenic_End":"end", "Accession_Number":"feature"})
    ori_df["type"] = "ori"
    ORI_OUTPUT = "ORI.summary.tsv"
    ori_df.to_csv(ORI_OUTPUT, sep="\t", index=False)
    count(ins_seq_df, ori_df, OUTPUT_NUM)


def singularize(OUTPUT_FILE):
    df = pd.read_csv("../filtered.ALL.plasmid_list.download.tsv", sep="\t")
    # assign single values for id and db for original plasmidscope db tsv file
    for i in range(0,len(df)):
        df.at[i, "Plasmid_ID"] = df.Plasmid_ID.iat[i].split(",")[0]
        df.at[i, "Data_Source"] = df.Data_Source.iat[i].split(",")[0]
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)

def final_df():
    # create a final df with all necessary info
    df = pd.read_csv("../final_df.tsv", sep="\t", index_col=False)
    df["AMR_binary"] = 0
    AMR_df = pd.DataFrame(columns=['Protein id', 'Contig id', 'Start', 'Stop', 'Strand', 'Element symbol', 'Element name', 'Scope', 'Type', 'Subtype', 'Class', 'Subclass', 'Method', 'Target length', 'Reference sequence length', '% Coverage of reference', '% Identity to reference', 'Alignment length', 'Closest reference accession', 'Closest reference name', 'HMM accession', 'HMM description'])
    PLASMIDS = df["Plasmid_ID"]

    # append AMRs
    for plasmid in PLASMIDS:
        path = f"../analysis/outputs/AMRs/{plasmid}.AMRs.tsv"
        if os.path.exists(path):
            amr = pd.read_csv(path, sep="\t")
            if not amr.empty:
                df.loc[df["Plasmid_ID"]==plasmid, "AMR_binary"] = 1
                AMR_df = pd.concat([AMR_df,amr], ignore_index=True)
    AMR_df.to_csv("../AMRs_df.tsv", sep="\t")
        
    df.to_csv("../final_df_AMRs.tsv", sep="\t", index=False)

def num_ins_seqs(FILE):
    is_df = pd.read_csv("IS.summary.tsv", sep="\t", index_col=False)
    main_df = pd.read_csv(FILE, sep="\t")
    main_df["num_IS"] = 0
    for id in is_df["seqID"].unique():
        sub_df = is_df[is_df["seqID"]==id]
        insseq = len(sub_df)
        main_df.loc[main_df["Plasmid_ID"]==id, "num_IS"] = insseq
    main_df.to_csv(FILE, sep="\t", index=False)

    




if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
