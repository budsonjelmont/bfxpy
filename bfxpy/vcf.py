import numpy as np
import pandas as pd
from collections import OrderedDict
from datetime import datetime as dt

from io import StringIO
import re

# Needed for df_to_pyvcf functionality:
#import vcf
#from vcf import utils

# Explode INFO columns into multiple columns named, e.g. 'INFODP', 'INFOFS', etc.
def explode_info_col(df):
  info_srs = df['INFO'].apply(lambda x: pd.DataFrame([{'INFO' + y.split('=')[0]:y.split('=')[1] for y in x.split(';')}]))
  info_df = pd.concat([r for r in info_srs]).reset_index(drop=True)
  df = df.join(info_df, how='left', lsuffix='', rsuffix='_info')
  df.drop(columns='INFO', inplace=True)
  return df

# Explode FORMAT columns into multiple columns named, e.g. 'FORMATDP', 'FORMATGQ', etc.
def explode_fmt_col(df, *sample_data_cols):
  for s in sample_data_cols:  # Currently only single sample VCFs are supported
    fmt_srs = df.apply(lambda x: pd.DataFrame([{'FORMAT' + y[0]:y[1] for y in zip(x['FORMAT'].split(':'),x[s].split(':'))}]),axis=1)
    fmt_df = pd.concat([r for r in fmt_srs]).reset_index(drop=True)
    df = df.join(fmt_df, how='left', lsuffix='', rsuffix='_format')
  df.drop(columns=['FORMAT'] + [s for s in sample_data_cols], inplace=True)
  return df

# Read VCF file and return a dataframe
def vcf_to_df(vcffile, *sample_data_cols, explode_fmt=False, explode_info=False):
  data = []
  #df = pd.read_csv(vcffile, sep='\t', comment='##', converters={'#CHROM':'CHROM'}) # pandas doesn't allow for comment indications
  with open(vcffile,'r') as f:
     # Change header line '##' to new comment character
     f1 = StringIO(re.sub(r'^##', '~', f.read(), flags=re.M))
  df = pd.read_csv(f1, sep='\t', comment='~', dtype={'#CHROM': str, 'POS': 'Int64', 'ID': str, 'REF': str, 'ALT': str})
  nan_columns = df.columns[pd.isna(df.columns)]
  df.rename(columns={np.nan:'dummyval'},inplace=True)
  df.columns = df.columns.fillna('SAMPLE_DATA')
  df.rename(columns={'#CHROM':'CHROM'}, inplace=True)
  cols=pd.Series(df.columns)
  for dup in cols[cols.duplicated()].unique(): 
    cols[cols[cols == dup].index.values.tolist()] = [dup + '.' + str(i) if i != 0 else dup for i in range(sum(cols == dup))]
  df.columns = cols
  if explode_info:
     df = explode_info_col(df)
#     info_srs = df['INFO'].apply(lambda x: pd.DataFrame([{'INFO' + y.split('=')[0]:y.split('=')[1] for y in x.split(';')}]))
#     info_df = pd.concat([r for r in info_srs]).reset_index(drop=True)
#     df = df.join(info_df, how='left', lsuffix='', rsuffix='_info')
#     df.drop(columns='INFO', inplace=True)
  if explode_fmt:
     df = explode_fmt_col(df, *sample_data_cols)
#    for s in sample_data_cols:  # Currently only single sample VCFs are supported
#      fmt_srs = df.apply(lambda x: pd.DataFrame([{'FORMAT' + y[0]:y[1] for y in zip(x['FORMAT'].split(':'),x[s].split(':'))}]),axis=1)
#      fmt_df = pd.concat([r for r in fmt_srs]).reset_index(drop=True)
#      df = df.join(fmt_df, how='left', lsuffix='', rsuffix='_format')
#    df.drop(columns=['FORMAT'] + sample_data_cols, inplace=True)
  return df

# Wrangle pandas dataframe into a consistent format
def preprocess_df(df,samplename):
    df.columns = df.columns.str.upper().str.replace('([FORMAT|INFO])_','\\1',regex=True)
    df.sort_values(['CHROM','POS','ID'],inplace=True)
    # Standardize data type
    df['CHROM'] = df['CHROM'].astype(str)
    df['POS'] = df['POS'].astype(int)
    # Gather up INFO and FORMAT columns
    if 'INFO' not in df.columns:
      infocols = list(df.columns[((df.columns.str.startswith('INFO')) & (df.columns!='INFO'))])
      df['INFO'] = df.apply(lambda x: ';'.join([y[4:] + '=' + str(x[y]) for y in infocols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)

    if 'FORMAT' not in df.columns or samplename not in df.columns:
      fmtcols = list(df.columns[((df.columns.str.startswith('FORMAT')) & (df.columns!='FORMAT'))])
      if 'FORMAT' not in df.columns:
        df['FORMAT'] = df.apply(lambda x: ':'.join([y[6:] for y in fmtcols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)
      if samplename not in df.columns:
        df[samplename] = df.apply(lambda x: ':'.join([str(x[y]) for y in fmtcols if ~pd.isna(x[y]) and x[y] is not None]),axis=1)

    # Add missing columns & fill nulls with default values
    optional_cols = {'FILTER':'PASS','QUAL':'.','INFO':'.','FORMAT':'.',samplename:'.'}
    missing_cols = optional_cols.keys() - set(df.columns)
    for col in missing_cols:
        df[col] = optional_cols[col]
    
    df.fillna(optional_cols,inplace=True)
    # Combine ALT1 & 2 into single field if there isn't already an ALT column defined
    if 'ALT' not in df.columns:
      df['ALT'] = df.apply(lambda x: ','.join([y for y in x[['ALT1','ALT2']] if ~pd.isna(y) and y is not None]),axis=1) 
    return df

# Chromosome sizes file generated using UCSC tools here https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# Example command:
# fetchChromSizes hg19 > hg19.chrom.sizes
def df_to_vcf(df, chr_size_file, output_file, build='GRCh38', source='df_to_vcf', samplename='SAMPLE_DATA'):
    df = preprocess_df(df,samplename)

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
        f.write('#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', samplename]) + '\n')

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

def trim_common_prefix(*sequences):
  if not sequences:
    return []
  reverses = [seq[::-1] for seq in sequences]
  rev_min = min(reverses)
  rev_max = max(reverses)
  if len(rev_min) < 2:
    return sequences
  for i, c in enumerate(rev_min[:-1]):
    if c != rev_max[i]:
      if i == 0:
        return sequences
      return [seq[:-i] for seq in sequences]
    return [seq[:-(i + 1)] for seq in sequences]

def df_to_pyvcf(df, l_trim=True, r_trim=True):
    # Preprocessing
    df.columns = df.columns.str.upper().str.replace('([FORMAT|INFO])_','\\1',regex=True)
    df.sort_values(['CHROM','POS','ID'],inplace=True)
    vcfreader = vcf.parser.Reader('nofile')
    reqdcols=['CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    missingcols = [col for col in reqdcols if col not in list(df.columns)]
    df[missingcols] = '.'
    # Troubleshooting: print all methods of vcf reader object:
    #object_methods = [method_name for method_name in dir(vcfreader) if callable(getattr(vcfreader, method_name))]
    #print(object_methods)
    # Special handling for qual column
    try:
      df['QUAL'] = df['QUAL'].astype(int)
    except ValueError as e:
      try:
        df['QUAL'] = df['QUAL'].astype(float)
      except ValueError as e:
        df['QUAL'] = None
    # Special handling for filter column
    # TODO enable handling for non-empty FORMAT column
    if (df['FORMAT'] == '.').all():
      df['FORMAT'] = None
    elif (df['FORMAT'] == '.').any():
      print('ERROR: Inconsistent entries in \'FORMAT\' field.')
    # Rewrite '.' in ALT field
    df['ALT'] = df['ALT'].replace({'.':None}).fillna(df['REF'])
    # Perform left- and right-trimmming, if requested
    print(df)
    if r_trim:
      df['REF'],df['ALT'] = zip(*df.apply(lambda x: utils.trim_common_suffix(x['REF'],x['ALT']), axis=1))
    if l_trim:
      rng = np.random.default_rng()
      ref_tmp_col = 'REF_' + str(rng.integers(low=0, high=77777777777, size=None, dtype=np.int64))
      alt_tmp_col = 'ALT_' + str(rng.integers(low=0, high=77777777777, size=None, dtype=np.int64))
      df[ref_tmp_col], df[alt_tmp_col] = zip(*df.apply(lambda x: utils.trim_common_suffix(x['REF'][::-1],x['ALT'][::-1]), axis=1))
      # Adjust POS for any rows that were left-trimmed
      df['POS'] = df['POS'] + (df['REF'].str.len() - df[ref_tmp_col].str.len())
      # Assign new REF & ALT alleles
      # TODO REF & ALT trimming below works, but vectorized doesn't
      df['REF'] = df.apply(lambda x: x['REF'][len(x['REF']) - len(x[ref_tmp_col]):], axis=1)
      df['ALT'] = df.apply(lambda x: x['ALT'][len(x['ALT']) - len(x[alt_tmp_col]):], axis=1)
      #df['REF'] = df['REF'].str[(df['REF'].str.len() - df[ref_tmp_col].str.len()):]
      #df['ALT'] = df['ALT'].str[(df['ALT'].str.len() - df[alt_tmp_col].str.len()):]
      # Re-sort dataframe
      df.sort_values('POS', ascending=True, inplace=True)
      print(df)
    return df.apply(lambda x: vcf.model._Record(x['CHR'], x['POS'], x['id'], x['REF'], [vcfreader._parse_alt(x['ALT'])], x['QUAL'], x['FILTER'], vcfreader._parse_info(x['INFO']), None, sample_indexes=None, samples=None), axis=1).values.tolist()

def df_to_jsonvcf(df, explode_info=True, explode_fmt=True, strip_non_alphanum=True):
  return df.to_dict(orient='records')