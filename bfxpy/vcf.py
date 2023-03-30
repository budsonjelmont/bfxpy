import pandas as pd
from collections import OrderedDict
from datetime import datetime as dt

# Chromosome sizes file generated using UCSC tools here https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# Example command:
# fetchChromSizes hg19 > hg19.chrom.sizes
def df_to_vcf(df, chr_size_file, output_file, build='GRCh38', source='df_to_vcf', samplename='FORMAT_DATA'):
    # Preprocessing
    df.columns = df.columns.str.upper().str.replace('([FORMAT|INFO])_','\\1',regex=True)
    df.sort_values(['CHROM','POS','ID'],inplace=True)

    # Fix data type issues
    df['CHROM'] = df['CHROM'].astype(str)

    # Gather up INFO and FORMAT columns
    infocols = list(df.columns[(df.columns.str.startswith('INFO'))])
    df['INFO'] = df.apply(lambda x: ';'.join([y[4:] + '=' + str(x[y]) for y in infocols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)

    fmtcols = list(df.columns[(df.columns.str.startswith('FORMAT'))])
    df['FORMAT'] = df.apply(lambda x: ':'.join([y[6:] for y in fmtcols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)
    df[samplename] = df.apply(lambda x: ':'.join([str(x[y]) for y in fmtcols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)

    # Add missing columns & fill nulls with default values
    optional_cols = {'FILTER':'PASS','QUAL':'.','INFO':'.','FORMAT':'.',samplename:'.'}
    missing_cols = optional_cols.keys() - set(df.columns)
    for col in missing_cols:
        df[col] = optional_cols[col]
    
    df.fillna(optional_cols,inplace=True)

    # Combine ALT1 & 2 into single field
    df['ALT'] = df.apply(lambda x: ','.join([y for y in x[['ALT1','ALT2']] if ~pd.isna(y) and y is not None]),axis=1)

    header = OrderedDict()
    header['fileformat'] = 'VCFv4.2'
    header['fileDate'] = dt.now().strftime('%Y%m%d')
    header['source'] = source
    header['reference'] = build
    header['contig'] = OrderedDict()

    # Get unique chromosome names and contig sizes
    chromosomes = df['CHROM'].unique()
    contig_sizes = {}
    for chrom in chromosomes:
        with open(chr_size_file, 'r') as f:
            for line in f:
                if line.startswith(chrom):
                    contig_sizes[chrom] = int(line.strip().split()[1])
                    break

        header['contig'][chrom] = OrderedDict()
        header['contig'][chrom]['ID'] = chrom
        header['contig'][chrom]['length'] = contig_sizes[chrom]
        header['contig'][chrom]['assembly'] = build
        #header['contig'][chrom]['md5'] = hashlib.md5(fasta.encode('utf-8')).hexdigest() # Not implemented, but maybe later

    # Write the header to the VCF file
    with open(output_file, 'w+') as f:
        for key, value in header.items():
            if key == 'contig':
                for subkey, subvalue in value.items():
                    contig_line = '##contig=<'
                    for k, v in subvalue.items():
                        contig_line += '{}={},'.format(k, v)
                    contig_line = contig_line.rstrip(',') + '>'
                    f.write(contig_line + '\n')
            else:
                f.write('##{}={}\n'.format(key, value))
        f.write('#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']) + '\n')

        # Write the variants to the VCF file
        for index, row in df.iterrows():
            chrom = row['CHROM']
            pos = row['POS']
            id = row['ID']
            ref = row['REF']
            alt = row['ALT']
            qual = row['QUAL']
            filt = row['FILTER']
            info = row['INFO']
            form = row['FORMAT']
            sampledata = row[samplename]

            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, id, ref, alt, qual, filt, info, form, sampledata))
