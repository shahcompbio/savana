import sys
import pandas as pd

if len(sys.argv) != 4:
    raise ValueError(f'\nUsage: python remap_patient_id.py <in.tsv> <in.metadata.paths.tsv> <out.tsv>\n')
in_tsv = sys.argv[1]
in_paths = sys.argv[2]
out_tsv = sys.argv[3]

pmap = pd.read_table(in_tsv) # isabl_patient_id -> patient
df = pd.read_table(in_paths)
id2patient = dict(zip(pmap['isabl_patient_id'], pmap['patient']))
df['isabl_patient_id'] = df['isabl_patient_id'].replace(id2patient)
df.to_csv(out_tsv, sep='\t', index=False)
