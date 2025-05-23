import numpy as np 
import pandas as pd 
import re 
import sys

# Arguments, make sure that the RGI table is the first argument, and the
# mobsuite database is second
rgi_path = sys.argv[1]
mob_path = sys.argv[2]

# Read both tables into memory
rgi_df = pd.read_table(rgi_path)
mob_df = pd.read_table(mob_path)

# RGI ouputs a string at the end of the contig IDs, e.g., _423. Thus, a regex
# is needed to convert the contig name back to its original, in order to
# properly join the results together.
new_contig = list()
for x in rgi_df['Contig']:
    y = re.sub(pattern='_[0-9]+\s{0,}$',
               repl='',
               string=x)
    new_contig.append(y)
rgi_df['Contig'] = new_contig


# Mob-suite reads the contig name the full string, but RGI ignores the string
# after the first whitespace. We will adjust this in the mob-suite output.
mob_new_contig = list()
for x in mob_df['contig_id']: 
    l = x.split()
    mob_new_contig.append(l[0])
mob_df['contig_id'] = mob_new_contig

# Pandas needs a common column name to join under, so we will change the name
# of the field in the mob results. Rename sample_id as well, it will move to
# the front of the ouput.
mob_df.rename(columns={'contig_id': 'Contig', 'sample_id': 'Sample'}, inplace=True)

# Before merging, make sure that the mobsuite contig names match up to the RGI
# contig names 
#if not set(rgi_df['Contig']).issubset(set(mob_df['Contig'])):
#    raise Exception(
"""
The contig names of the MOB-suite results are not a subset of the contig 
names of the RGI results.
"""
#)

# Merge the results by performing a left join. RGI results appear first in the
# final table. There should be no duplication of Contigs in the mobDB, but
# there may be many ORFs per contig in the RGI results, so we set the
# many_to_one validation method.
merged_df = pd.merge(left=rgi_df, 
                     right=mob_df, 
                     on='Contig', 
                     how='left', 
                     validate='many_to_one')

# Add the Sample column from the mob-suite output to the very front, it just
# feels more natural.
col = merged_df.pop("Sample")
merged_df.insert(0, col.name, col)

# Output a file
merged_df.to_csv(path_or_buf='merged_tables.csv', 
                 sep=',', 
                 index=False)

