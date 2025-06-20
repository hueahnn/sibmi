import pandas as pd
t = pd.read_tsv("ALL.plasmid_list.download.tsv")
filtered = t[(t["Topology"] == "circular") & (t["Completeness"] == "complete")]
PLASMIDS = filtered["Plasmid_ID"]
rule all:
    input: expand("results/{plasmid}/{plasmid}.gbff", plasmid = PLASMIDS)
rule annotate:
    input: 
        "fasta-inputs-2/{plasmid}"
    output:
        "results/{plasmid}/{junk}.gbff"
    conda: "bakta.yaml"
    shell:
        "bakta --db /n/data1/hms/dbmi/baym/databases/bakta_dbv6/db {input} --output {output} --meta --force"
