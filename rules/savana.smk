import os
import glob
import pysam

def _get_system_id_from_bam(bam_path):
    #bam_path = wildcards.input.tumor_bam
    bam = pysam.AlignmentFile(bam_path)
    bam_header_dict = bam.header.as_dict()
    header_last_CL = bam_header_dict['PG'][-1]['CL'].split(' ')[-1]
    system_id = header_last_CL.replace('.bam', '')
    return system_id

def _get_ont_normal_bam(wildcards):
    if wildcards.sample.startswith('PMC-01'):
        return '/juno/work/shah/isabl_data_lake/analyses/84/23/38423/results/minimap2/SHAH_H003452_T01_01_WG01_R1.sorted.bam' # PMN-01-N1
    #paths_path = "/juno/work/shah/users/chois7/tables/paths.ONT-NANOSEQ.tsv"
    paths_path = config["paths_file"]
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_patient_id']==wildcards.sample]
    path = paths[(paths["sample_category"] == "NORMAL") & (paths["result_type"] == "sorted_bam")].reset_index()["result_filepath"][0]
    assert os.path.exists(path), path
    return path

def _get_ont_tumor_bam(wildcards):
    paths_path = config["paths_file"]
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_patient_id']==wildcards.sample]
    path = paths[(paths["sample_category"] == "TUMOR") & (paths["result_type"] == "sorted_bam")].reset_index()["result_filepath"][0]
    assert os.path.exists(path), path
    return path

out_dir = config["out_dir"]

rule run_savana:
    input: 
        normal_bam = _get_ont_normal_bam,
        tumor_bam = _get_ont_tumor_bam,
    output:
        #vcf = lambda wildcards: f'results/{{sample}}/savana/{_get_system_id_from_bam(wildcards.input)}.',
        #vcf = "results/{sample}/savana/variant.stats",
        vcf= os.path.join(out_dir, "{sample}/savana/variant.stats"),
    #singularity: "docker://soymintc/savana:latest"
    params:
        #outdir = "results/{sample}/savana",
        outdir = os.path.join(out_dir,"{sample}/savana"),
        ref = config['ref'],
        ref_index = config['ref_index'],
    threads: 12
    shell:
        """
        savana --tumour {input.tumor_bam} --normal {input.normal_bam} --ref {params.ref} --ref_index {params.ref_index} --outdir {params.outdir} --threads {threads} --mapq 0
        """

def _get_savana_lenient_vcf(wildcards):
    vcfs = glob.glob(f'results/{wildcards.sample}/savana/*.somatic.sv_breakpoints.lenient.vcf')
    assert len(vcfs) == 1, vcfs
    vcf = vcfs[0]
    return vcf

rule conform_savana_svs:
    input:
        vcf = _get_savana_lenient_vcf,
        stats = "results/{sample}/savana/variant.stats",
    output:
        tsv = "results/{sample}/svs/{sample}.svs.tsv",
    shell:
        """
        python scripts/conform_savana_svs.py -i {input.vcf} -o {output.tsv}
        """

rule annotate_genes:
    input:
        tsv = "results/{sample}/svs/{sample}.svs.tsv",
    output:
        tsv = "results/{sample}/svs/{sample}.svs.genes.tsv",
    params:
        fai = config['genome_fai'],
        gtf = config['genome_gtf_gz'],
    shell:
        """
        python scripts/annotate_genes.py -s {input.tsv} -o {output.tsv} --gtf {params.gtf} --fai {params.fai}
        """

