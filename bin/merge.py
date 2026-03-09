import pandas as pd
import re
import sys

# Improved contig normalization
def normalize_contig(name):
    # Take first whitespace-separated part
    name = name.split()[0]
    # Remove coverage info
    name = re.sub(r'_cov_[0-9.]+', '', name)
    # Remove multiple pilon suffixes
    name = re.sub(r'(_pilon)+', '', name)
    # Keep NODE_X_length_YYYYYY intact, remove trailing numbers after length
    # Match pattern: NODE_<num>_length_<num>_optionalExtra -> NODE_<num>_length_<num>
    m = re.match(r'(NODE_\d+_length_\d+)', name)
    if m:
        name = m.group(1)
    return name

# Arguments: RGI table, MOB table, Sample name
rgi_path = sys.argv[1]
mob_path = sys.argv[2]
sample = sys.argv[3]



# Read tables
rgi_df = pd.read_table(rgi_path)
mob_df = pd.read_table(mob_path)


rgi_df['Contig'] = rgi_df['Contig'].apply(normalize_contig)
mob_df['contig_id'] = mob_df['contig_id'].apply(normalize_contig)
mob_df.rename(columns={'contig_id': 'Contig'}, inplace=True)
#mob_df['contig_id'] = mob_df['contig_id'].apply(normalize_contig)
#print("RGI contigs:", set(rgi_df['Contig']))
#print("MOB contigs:", set(mob_df['Contig']))
#print("Intersection:", set(rgi_df['Contig']) & set(mob_df['Contig']))
#print(set(mob_df['Contig']) & set(rgi_df['Contig']))
# Rename MOB contig column for merge
#mob_df.rename(columns={'contig_id': 'Contig'}, inplace=True)

# Merge RGI and MOB: keeps all MOB columns
merged_df = pd.merge(
    left=rgi_df,
    right=mob_df,
    on='Contig',
    how='left',         # keeps all RGI rows
    validate='many_to_one'
)

# Insert sample name
merged_df.insert(0, 'Sample', sample)

# Save output
merged_df.to_csv('merged_tables.csv', sep=',', index=False)
