import os
PLASMIDS = os.listdir("fasta-inputs")
rule all:
    input: expand("results/{plasmid}", plasmid = PLASMIDS)
rule annotate:
    input: 
        "fasta-inputs-2/{plasmid}"
    output:
        "results/{plasmid}/{plasmid}.gbff"
    conda: "bakta.yaml"
    shell:
        "bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db {input} --output $(dirname {output}) --meta --force"
