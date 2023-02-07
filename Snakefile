import os
import pandas as pd

configfile: "config.yaml"
patientlist = pd.read_csv(config["patientlist"])
SAMPLES = patientlist["isabl_patient_id"]

log_dir = config['log_dir']
tmp_dir = config['tmp_dir']
for dir_path in [log_dir, tmp_dir]:
    if not os.path.exists(dir_path):
        os.system(f"mkdir {dir_path}")

rule all:
    input: 
        expand("results/{sample}/savana/somatic.sv_breakpoints.strict.vcf", sample=SAMPLES),

include: "rules/savana.smk"
