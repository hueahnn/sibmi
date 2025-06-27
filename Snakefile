import os
PATH = "first-week/fasta-inputs-2"
PLASMIDS = [f.replace(".fasta", "") for f in os.listdir(PATH) if f.endswith(".fasta")]

def main_dir(child):
    return os.path.join(PATH, child)

rule all:
    input: expand("first-week/results/{plasmid}/{plasmid}.gbff", plasmid = PLASMIDS)
rule annotate:
    input: 
        "first-week/fasta-inputs-2/{plasmid}.fasta"
    output:
        "first-week/results/{plasmid}/{plasmid}.gbff"	
    conda: "conda-envs/bakta.yaml"
    shell:
        "bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db {input} --output $(dirname {output}) --meta --force"
