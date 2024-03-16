#!/usr/bin/env python


import pandas as pd
import sys
import csv
import numpy as np

def convert_gtf_to_bed12(gtf_path):
    # Define column names for the GTF file
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    
    # Load the GTF file
    gtf_df = pd.read_csv(gtf_path, sep='\t', comment='#', names=gtf_columns, low_memory=False)
    
    # Filter for exons and CDS
    exon_df = gtf_df[gtf_df['feature'].isin(['exon', 'CDS'])].copy()
    
    # Extract gene_id and transcript_id from attributes
    exon_df['gene_id'] = exon_df['attributes'].str.extract('gene_id "([^"]+)"')
    exon_df['transcript_id'] = exon_df['attributes'].str.extract('transcript_id "([^"]+)"')
    
    # Determine thickStart and thickEnd from CDS regions
    cds_df = exon_df[exon_df['feature'] == 'CDS'].copy()
    thick_df = cds_df.groupby('transcript_id')[['start', 'end']].agg({'start': 'min', 'end': 'max'}).rename(columns={'start': 'thickStart', 'end': 'thickEnd'})

    # Focus on exons for the rest
    exon_df = exon_df[exon_df['feature'] == 'exon']
    
    # Merge thickStart and thickEnd information back into exon_df
    exon_df = exon_df.merge(thick_df, on='transcript_id', how='left')
    
    # Now, use the start of the first exon as chromStart and the end of the last exon as chromEnd
    # thickStart and thickEnd are based on CDS, but if there's no CDS (NaN), use exon positions
    bed12_df = exon_df.groupby('transcript_id').apply(lambda x: pd.Series({
        'chrom': x['seqname'].iloc[0],
        'chromStart': x['start'].min() - 1,
        'chromEnd': x['end'].max(),
        'name': x['transcript_id'].iloc[0],
        'score': 0,
        'strand': x['strand'].iloc[0],
        'thickStart': np.nanmin(x['thickStart'].fillna(x['start']).values) - 1,  # Use CDS start if available, else exon start
        'thickEnd': np.nanmax(x['thickEnd'].fillna(x['end']).values),  # Use CDS end if available, else exon end
        'itemRgb': '0,0,0',
        'blockCount': len(x),
        'blockSizes': ','.join((x['end'] - x['start'] + 1).astype(str)),
        'blockStarts': ','.join((x['start'] - x['start'].min()).astype(str)),
    })).reset_index(drop=True)
    
    # Select and reorder columns to match BED12 format
    bed12_columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                     'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    bed12_df = bed12_df[bed12_columns]
    
    bed12_df.sort_values(by=['chrom', 'chromStart', 'chromEnd', 'strand', 'name'], inplace=True)
    
    return bed12_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <path_to_gtf_file>")
        sys.exit(1)
    
    gtf_file_path = sys.argv[1]
    bed12_df = convert_gtf_to_bed12(gtf_file_path)
    
    # Example to save the output, adjust path as necessary
    bed12_df.to_csv(gtf_file_path.replace('.gtf', '.bed'), sep='\t', index=False, header=False)

