import re
import argparse
import pandas as pd
from collections import defaultdict
pd.options.mode.chained_assignment = None  # default='warn'; this is to try to suppress a useless warning

def parse_args():
    parser = argparse.ArgumentParser(description='Process SAVANA vcf')
    parser.add_argument('-i', '--input', help='Input SAVANA vcf file', required=True)
    parser.add_argument('-o', '--output', help='Output adjacency tsv file', required=True)
    # parser.add_argument('-s', '--sample', help='Input sample ID for vcf header', required=True)
    args = parser.parse_args()
    return args

class Breakpoint:
    def __init__(self, brk):
        self.chrom = brk['chrom']
        self.pos = brk['pos']
        self.self_id = brk['breakend'] # ID
        brk_ix = int(self.self_id[-1])
        assert brk_ix in {1, 2}, f'ERROR: brk_ix = {brk_ix}'
        self.ref = brk['ref']
        self.alt = brk['alt']
        self.info = brk['INFO']
        bp_notation = re.search('BP_NOTATION=([^;]+);', self.info).groups()[0]
        self.adjacency_id = brk['adjacency'] # ID
        self.mate_id = f'{self.adjacency_id}_{3-brk_ix}'
        
        _strand_combination = {'++', '--', '+-', '-+'}
        if bp_notation == '<INS>':
            self.strand = None
        elif bp_notation in _strand_combination:
            self.strand = bp_notation[brk_ix-1]
        else:
            raise ValueError(f'ERROR: bp_notation = {bp_notation}')

class Adjacency:
    def __init__(self, brks): # brks <- paired dataframe
        assert brks.shape[0] == 2, brks
        brk1, brk2 = brks.iloc[0], brks.iloc[1]
        self.brk1 = Breakpoint(brk1)
        self.brk2 = Breakpoint(brk2)
        self.type = 'n/a'
        
        self.type = self.get_svtype()
        self.length = abs(self.brk2.pos - self.brk1.pos)
        if self.type == 'translocation': 
            self.length = int(3e9)
    
    def get_svtype(self): # N-> <-N // <-N N-> // N<- N<- // ->N ->N
        if self.brk1.chrom != self.brk2.chrom: # - <TRA>
            return 'translocation'
        if (self.brk1.strand, self.brk2.strand) == ('+', '+'):
            return 'inversion'
        elif (self.brk1.strand, self.brk2.strand) == ('-', '-'):
            return 'inversion'
        elif (self.brk1.strand, self.brk2.strand) == ('+', '-'):
            return 'deletion'
        elif (self.brk1.strand, self.brk2.strand) == ('-', '+'):
            return 'duplication'
        else:
            raise ValueError(f'ERROR: (strand1, strand2) = ({self.brk1.strand}, {self.brk1.strand})')

def resolve_breakpoints(df):
    df = df.copy()
    chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y']
    svsv_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2',
                 'type', 'length']
    svsv = pd.DataFrame(columns=svsv_cols) # savana sv

    df['adjacency'] = df['breakend'].str.rsplit('_', 1, expand=True)[0]
    
    svtype_cnt = defaultdict(int)
    for i, (_, brks) in enumerate(df.groupby('adjacency')):
        brks_in_adj = brks.shape[0]
        if brks_in_adj == 2:
            brk1, brk2 = brks.iloc[0], brks.iloc[1]
            brk1.chrom = brk1.chrom.replace('chr', '')
            brk2.chrom = brk2.chrom.replace('chr', '')
            if brk1.chrom not in chroms: continue
            if brk2.chrom not in chroms: continue
            if brk1.chrom == brk2.chrom:
                assert brk1.pos <= brk2.pos, (brk1.pos, brk2.pos)
            adj = Adjacency(brks)
            svtype_cnt[adj.type] += 1
            line = [adj.brk1.chrom, adj.brk1.pos, adj.brk1.strand, adj.brk2.chrom, adj.brk2.pos, adj.brk2.strand, adj.type, adj.length]
        elif brks_in_adj == 1:
            brk = brks.squeeze()
            brk1 = Breakpoint(brk)
            if brk1.chrom not in chroms: continue
            assert brk['alt'] == '<INS>', brk
            svtype = 'insertion'
            match = re.search('INSSEQ=([A-Z]+);', brk['INFO'])
            insseq = match.groups()[0]
            svlength = len(insseq)
            svtype_cnt[svtype] += 1
            line = [brk1.chrom, brk1.pos, brk1.strand, brk1.chrom, brk1.pos, brk1.strand, svtype, svlength]
        svsv.loc[i] = line

        # N-> <-N // <-N N-> // N-> N-> // <-N <-N
    # svsv['type'] = svsv['type'].replace('insertion', 'ins').replace('deletion', 'del').replace('inversion', 'inv').replace('duplication', 'dup')
    return svsv 

def get_svs(svs_path, sv_types={'inversion', 'translocation'}):
    vcf_info = []
    for line in open(svs_path, 'r'):
        if line.startswith('##'):
            vcf_info.append(line)

    vcf_cols = ['chrom', 'pos', 'breakend', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample']

    svs = pd.read_table(svs_path, comment='#', names=vcf_cols)
    svs = resolve_breakpoints(svs)
    # svs = svs[
    #     svs['type'].isin(sv_types) # filter in INV TRA only
    # ]
    return vcf_info, svs


if __name__ == "__main__":
    args = parse_args()
    vcf_info, svs = get_svs(args.input)
    svs.to_csv(args.output, sep='\t', index=False)
    
