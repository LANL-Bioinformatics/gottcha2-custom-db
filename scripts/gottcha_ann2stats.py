#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os

asr=""

try:
    asr=sys.argv[1]
except:
    print( "\nUSAGE: %s [bwa-index].ann\n"%sys.argv[0] )
    sys.exit(1)

sys.stderr.write( "[INFO] Loading %s..."%asr )
#chunksize = 10**6 
df = pd.read_csv(
    asr,
    sep='|',
    header=None,
    names=['name','start','end','taxid','desc'],
    usecols=[1,2,3],
    dtype={'taxid':np.str}
    #verbose=True
    #chunksize=chunksize
)
sys.stderr.write( "Done.\n" )

sys.stderr.write( "[INFO] Dropping empty rows..." )
df = df[df.taxid.notnull()]
sys.stderr.write( "Done.\n" )

sys.stderr.write( "[INFO] Aggregating..." )
df['len']=df.end-df.start+1
stats_df = df.groupby('taxid').agg({'len': ['min', 'max', 'sum', 'count']}).len
stats_df.reset_index(inplace=True)
sys.stderr.write( "Done.\n" )

stats_df['Rank']='strain'
stats_df['Name']=''
stats_df['Superkingdom']=''
stats_df['min']=stats_df['min'].astype('int64')
stats_df['max']=stats_df['max'].astype('int64')
stats_df['sum']=stats_df['sum'].astype('int64')

stats_df.rename(
    index=str,
    inplace=True,
    columns={"taxid": "Taxid", "min": "Min", "max": "Max", "count": "NumOfSeq", "sum": "TotalLength"}
)

sys.stderr.write( "[INFO] Output stats..." )
stats_df.to_csv(
    os.sys.stdout,
    sep="\t",
    header=True,
    index=False,
    columns=['Rank','Name','Taxid','Superkingdom','NumOfSeq','Max','Min','TotalLength']
)
sys.stderr.write( "Done.\n" )
