import sys
import argparse
import pybedtools
import tabix
import numpy as np
import pandas as pd

def _get_svs(svs_path):
    svs = pd.read_table(svs_path)
    return svs

def parse_args():
    parser = argparse.ArgumentParser(description='Annotate SVs as tsv table')
    parser.add_argument('-s', '--svs', help='Input SAVANA SV tsv', required=True)
    parser.add_argument('--gtf', help='Input reference gtf', required=True)
    parser.add_argument('--fai', help='Input reference genome fasta index', required=True)
    parser.add_argument('-o', '--out_svs', help='Output SVs with gene annotations', required=True)
    args = parser.parse_args()
    return args


def _extract_gene_name(gene_info):
    # gene_info[0]
    # 'gene_id "ENSG00000227160.3"; gene_type "transcribed_unitary_pseudogene"; gene_name "THEM7P"; level 2; hgnc_id "HGNC:50386"; havana_gene "OTTHUMG00000154120.3";'
    gene_names = [g.split('gene_name "')[1].split('";')[0] for g in gene_info]
    gene_names = list(set(gene_names))
    gix = 0
    if len(gene_names) > 1:
        print(f'{gene_info}: more than one gene_names:{gene_names}', file=sys.stderr)
        gene_names_clean = [g for g in gene_names if g.count('.') == 0 and g.count('-') == 0]
        if len(gene_names_clean) >= 1:
            gix = gene_names.index(gene_names_clean[0])
    elif len(gene_names) == 0:
        return ''
    gene_name = gene_names[gix]
    return gene_name

def get_insertion_genes(gtf, brk, fai_path):
    """Get genes near an insertion with strand taken into account

    :param gtf: BedTools GTF
    :type gtf: pybedtools.bedtool.BedTool
    :param brk: (chrom, pos, strand) tuple
    :type gtf: tuple
    :param fai_path: path to fai file
    :type fai_path: str
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    chrom, pos, strand = brk
    canon_chroms = set([f'chr{c}' for c in range(1,22+1)] + ['chrX', 'chrY'])
    if chrom not in canon_chroms: 
        return np.nan
    bed = pybedtools.BedTool(f'{chrom} {pos-1} {pos}', from_string=True)
    intersect = bed.intersect(gtf, g=fai_path, 
                              wb=True) # report 'b[gtf]' records
    gene_info = [g[-1] for g in intersect]
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def get_overlapping_genes(tbx, brk):
    """Get genes overlapping a breakpoint 

    :param tbx: PyTabix GTF
    :type tbx: pybedtools.bedtool.BedTool
    :param brk: (chrom, pos, strand) tuple
    :type brk: tuple
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    chrom, pos, _ = brk
    pos = int(pos)
    results = tbx.query(chrom, pos-1, pos)
    gene_info = [g[8] for g in results]
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def get_closest_genes(gtf, brk, fai_path):
    """Get genes near a breakpoint with strand taken into account

    :param gtf: BedTools GTF
    :type gtf: pybedtools.bedtool.BedTool
    :param brk: (chrom, pos, strand) tuple
    :type gtf: tuple
    :param fai_path: path to fai file
    :type fai_path: str
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    chrom, pos, strand = brk
    canon_chroms = set([f'chr{c}' for c in range(1,22+1)] + ['chrX', 'chrY'])
    if chrom not in canon_chroms: 
        return ''
    bed = pybedtools.BedTool(f'{chrom} {pos-1} {pos}', from_string=True)
    closest = bed.closest(gtf, g=fai_path,
                          fu=(strand=='+'), # get upstream if '+'
                          fd=(strand=='-'), # get downstream if '-'
                          D="a") # report distance from "a[bed]"=>brk        
    gene_info = [g[11] for g in closest]
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def add_gene_names_to_svs(svs, gtf, tbx, fai_path):
    gene_names = {1:[], 2:[]}
    svs = svs.copy()
    for rix, row in svs.iterrows():
        if rix % 10 == 0: print(f'progress: {rix}/{svs.shape[0]}')
        chrom1, pos1, strand1, chrom2, pos2, strand2, svtype, svlength = row.tolist()
        assert (svtype=='translocation') or (svtype!='translocation' and chrom1==chrom2 and pos1<=pos2), row
        brks = [(chrom1, pos1, strand1), (chrom2, pos2, strand2)]
        if svtype == 'insertion':
            brk = brks[0]
            gene_name = get_insertion_genes(gtf, brk, fai_path)
            gene_names[1].append(gene_name)
            gene_names[2].append(gene_name)
        else:
            for ix, brk in enumerate(brks):
                bix = ix + 1
                gene_name = get_overlapping_genes(tbx, brk)
                if len(gene_name) == 0: # ''
                    gene_name = get_closest_genes(gtf, brk, fai_path)
                gene_names[bix].append(gene_name)

    for ix in [1, 2]:
        svs[f'gene_name_{ix}'] = gene_names[ix]
    return svs

if __name__ == "__main__":
    args = parse_args()
    fai_path = args.fai
    gtf = pybedtools.BedTool(args.gtf)
    tbx = tabix.open(args.gtf)
    svs = _get_svs(args.svs)
    gsvs = add_gene_names_to_svs(svs, gtf, tbx, fai_path)
    gsvs.to_csv(args.out_svs, sep='\t', index=False)
