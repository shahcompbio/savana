def _get_ont_normal_bam(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.ONT-NANOSEQ.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_patient_id']==wildcards.sample]
    path = paths[(paths["sample_category"] == "NORMAL") & (paths["result_type"] == "sorted_bam")].reset_index()["result_filepath"][0]
    return path

def _get_ont_tumor_bam(wildcards):
    paths_path = "/juno/work/shah/users/chois7/tables/paths.ONT-NANOSEQ.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_patient_id']==wildcards.sample]
    path = paths[(paths["sample_category"] == "TUMOR") & (paths["result_type"] == "sorted_bam")].reset_index()["result_filepath"][0]
    return path

rule run_savana:
    input: 
        normal_bam = _get_ont_normal_bam,
        tumor_bam = _get_ont_tumor_bam,
    output:
        outdir = "results/{sample}/savana",
        done = "results/{sample}/savana/DONE",
    singularity: "docker://soymintc/savana:latest"
    params:
        ref = config['ref'],
        ref_index = config['ref_index'],
    threads: 4
    shell:
        """
        savana --tumour {input.tumor_bam} --normal {input.normal_bam} --ref {params.ref} --ref_index {params.ref_index} --outdir {output.outdir} --threads {threads}
        touch {output.done}
        """
