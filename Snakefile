import os
import pandas as pd

configfile: "config.yaml"
#patientlist = pd.read_csv(config["patientlist"])
#SAMPLES = patientlist["isabl_patient_id"]
#SAMPLES = ['PMC-01-PDOX1', 'PMC-01-PDOX2']
SAMPLES=["AK-RT-004"]

log_dir = config['log_dir']
tmp_dir = config['tmp_dir']
out_dir = config["out_dir"]
for dir_path in [log_dir, tmp_dir]:
    if not os.path.exists(dir_path):
        os.system(f"mkdir {dir_path}")

rule all:
    input: 
        expand(os.path.join(out_dir, "{sample}/savana/variant.stats"), sample=SAMPLES),
        #expand("results/{sample}/svs/{sample}.svs.genes.tsv", sample=SAMPLES), # for PMC-01

include: "rules/savana.smk"
