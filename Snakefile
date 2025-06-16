import os
PLASMIDS = os.listdir("fasta-inputs")
print(PLASMIDS)
rule all:
    input: expand("results/{plasmid}", plasmid = PLASMIDS)
rule annotate:
    input: 
        "fasta-inputs/{plasmid}"
    output:
        directory("results/{plasmid}")
    conda: "bakta.yaml"
    shell:
        "bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db {input} --output {output}"
