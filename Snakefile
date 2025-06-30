import os
PATH = "bakta-annotations"
INPUT_FILE = "filtered.plasmid.fasta.txt"
IDS = []
with open(INPUT_FILE, 'r') as f:
    IDS = [line.strip() for line in f if line.strip()]

def main_dir(child):
    return os.path.join(PATH, child)

rule all:
    input: expand("{path}/results/{plasmid}/{plasmid}.gbff", path = PATH, plasmid = IDS)
rule annotate:
    input:
        lambda wildcards: f"{wildcards.path}/../fasta-files/{wildcards.plasmid.split('_')[0]}/fasta/{wildcards.plasmid}.fasta"
    output:
        "{path}/results/{plasmid}/{plasmid}.gbff"
    conda: "conda-envs/bakta.yaml"
    shell:
        "bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db {input} --output $(dirname {output}) --meta --force"
