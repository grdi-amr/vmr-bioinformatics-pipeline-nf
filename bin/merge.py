import pandas as pd
import re
import sys

# Improved contig normalization: match by NODE ID only
def normalize_contig(name):
    name = name.split()[0]

    # Case 1: NODE format (SPAdes)
    m = re.match(r'(NODE_\d+)', name)
    if m:
        return m.group(1)

    # Case 2: contig format (Prokka / MOB-suite)
    m = re.match(r'(contig\d+)', name)
    if m:
        return m.group(1)

    return name

# Arguments: RGI table, MOB table, Sample name
rgi_path = sys.argv[1]
mob_path = sys.argv[2]
sample = sys.argv[3]

# Read tables
rgi_df = pd.read_table(rgi_path)
mob_df = pd.read_table(mob_path)

# Normalize contigs
rgi_df['Contig'] = rgi_df['Contig'].apply(normalize_contig)
mob_df['contig_id'] = mob_df['contig_id'].apply(normalize_contig)
mob_df.rename(columns={'contig_id': 'Contig'}, inplace=True)

print("RGI contigs:", set(rgi_df['Contig']))
print("MOB contigs:", set(mob_df['Contig']))
print("Intersection:", set(rgi_df['Contig']) & set(mob_df['Contig']))

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
