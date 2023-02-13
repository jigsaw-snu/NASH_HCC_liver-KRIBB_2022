import os
import pandas as pd
import gseapy as gp
from statsmodels.robust import mad


def CountExtractor(file_path: str) -> tuple:
    """
        Extract Count Values from RSEM results
    """

    cnt_df = pd.DataFrame()
    tpm_df = pd.DataFrame()
    fpkm_df = pd.DataFrame()

    # list of RSEM results
    rsem_results = os.listdir(file_path)

    first_pass = True
    for rsem_res in rsem_results:
        tmp_res = pd.read_csv(os.path.join(file_path, rsem_res), sep='\t')

        # remove version from gene ID
        tmp_res['gene_id'] = pd.Series(map(lambda x: x.split('.')[0], tmp_res['gene_id']))

        tmp_res.index = tmp_res['gene_id']

        if first_pass:
            # explicitly convert Pandas Series (Single Column) to Pandas DataFrame
            len_df = tmp_res['length'].to_frame()
            first_pass = False
        else:
            if all(tmp_res['index'] != cnt_df.index):
                raise Exception("ERR : Gene Index of current RSEM result is different from previous one")
        
        # CCI4_1x_6.genes.results --> CCI4_1x_6
        sample_name = rsem_res[:-14]

        cnt_df[sample_name] = tmp_res['expected_count']
        tpm_df[sample_name] = tmp_res['TPM']
        fpkm_df[sample_name] = tmp_res['FPKM']

    return (cnt_df, tpm_df, fpkm_df, len_df)



def MADFilter(df: pd.DataFrame) -> pd.DataFrame:
    """
        Remove Outliers using MAD criteria
        @ MAD criteria : 
            -> Remove Genes with value of (median - 3 * MAD)
    """

<<<<<<< Updated upstream
    upper_bound = df.median(axis=0) + 3*mad()
=======
    upper_bound = df.median(axis=0) + 3*mad(df, axis=0)  # median(axis=0) : median for each column vector
    lower_bound = df.median(axis=0) - 3*mad(df, axis=0)  # mad(axis=0) : MAD for each column vector
    
    # retain genes that satisfies MAD criteria for every sample
    return df[((df >= lower_bound) & (df < upper_bound)).sum(axis=1) == df.shape[1]]  # sum(axis=1) : sum for each row vector (count if it's True for every sample)
>>>>>>> Stashed changes



def CalculateTPM(cnt_df: pd.DataFrame, len_df: pd.DataFrame):
    """
        Calculate TPM values from unnormalized count data
    """
    
    common_idx = set(cnt_df.index) & set(len_df.index)
    cnt = cnt_df.loc[list(common_idx), :]
    length = len_df.loc[list(common_idx), :]

    rpk = cnt.div(length.values, axis=0) * 1e3
    
    tpm = rpk.div(rpk.sum().values) * 1e6
    
    return tpm



def FindHighVarGenes():
    """
        Find genes that are highly variable within each samples
        Two version of outputs exist:
            (1) 
            (2)
    """



def Ens2Sym(df: pd.DataFrame) -> pd.DataFrame:
    """
        Convert Ensembl Gene ID into Gene Symbols
    """

    # make sure that rownames are Ensembl Gene ID
    if 'ENS' not in df.index.tolist()[0]:
        raise Exception("ERR : rownames of given dataframe must be an Ensembl Gene ID")

    # since biomart only takes so tiny amount of queries (upto 300 queries), we need batch-work
    batch_size = 300
    nr_batches = (df.shape[0] // batch_size) + ((df.shape[0] % batch_size) > 0)

    # empty dataframe for final result. column names need to be same as biomart output for merging later
    conversion_res = pd.DataFrame(columns=['ensembl_gene_id', 'external_gene_name'])  

    # biomart API from gseapy
    bm = gp.Biomart()

    for i in range(nr_batches):
        queries = {'ensembl_gene_id' : df.index.tolist()[batch_size * i : batch_size * (i+1)]}

        conversion_batch_res = bm.query(dataset='mmusculus_gene_symbol',
                                        attributes=['ensembl_gene_id', 'external_gene_name'],
                                        filters=queries)
        
        # concatenate batch result to dataframe of final result
        conversion_res = pd.concat([conversion_res, conversion_batch_res])
    
    # remove rows with unconverted gene ID
    conversion_res = conversion_res.dropna(axis=0)

    # merge original dataframe and conversion result (make converted ID a new column of original dataframe)
    merged_res = df.reset_index(inplace=False)  # Two merging dataframes must have same index names
    merged_res = merged_res.rename(columns={'gene_id':'ensembl_gene_id'})  # Two merging dataframe must have same name for a column where mergence occurs
    merged_res = pd.merge(merged_res, conversion_res, how='inner', on='ensembl_gene_id')  # merge conversion result and original dataframe
    merged_res = merged_res.drop(columns=['ensembl_gene_id'])  # we do not need original ensembl gene ID anymore
    merged_res = merged_res.set_index('external_gene_name')  # set converted gene ID as a row index

    # handle multiple entries (in case where one ensembl ID mapped to multiple entrez ID)
    merged_res = merged_res.groupby(merged_res.index).mean().round()  # set averaged and rounded-down value for a duplicated rows
    merged_res.sort_index(axis=1, inplace=True)  # sort column names (sample names) alphabetically

    return merged_res



def Ens2Ent(df: pd.DataFrame) -> pd.DataFrame:
    """
        Convert Ensembl Gene ID into Entrez ID
    """

    # make sure that rownames are Ensembl Gene ID
    if 'ENS' not in df.index.tolist()[0]:
        raise Exception("ERR : rownames of given dataframe must be an Ensembl Gene ID")

    # since biomart only takes so tiny amount of queries (upto 300 queries), we need batch-work
    batch_size = 300
    nr_batches = (df.shape[0] // batch_size) + ((df.shape[0] % batch_size) > 0)

    # empty dataframe for final result. column names need to be same as biomart output for merging later
    conversion_res = pd.DataFrame(columns=['ensembl_gene_id', 'entrezgene_id'])  

    # biomart API from gseapy
    bm = gp.Biomart()

    for i in range(nr_batches):
        queries = {'ensembl_gene_id' : df.index.tolist()[batch_size * i : batch_size * (i+1)]}

        conversion_batch_res = bm.query(dataset='mmusculus_gene_symbol',
                                        attributes=['ensembl_gene_id', 'entrez_gene_id'],
                                        filters=queries)
        
        # concatenate batch result to dataframe of final result
        conversion_res = pd.concat([conversion_res, conversion_batch_res])
    
    # remove rows with unconverted gene ID
    conversion_res = conversion_res.dropna(axis=0)

    # merge original dataframe and conversion result (make converted ID a new column of original dataframe)
    merged_res = df.reset_index(inplace=False)  # Two merging dataframes must have same index names
    merged_res = merged_res.rename(columns={'gene_id':'ensembl_gene_id'})  # Two merging dataframe must have same name for a column where mergence occurs
    merged_res = pd.merge(merged_res, conversion_res, how='inner', on='ensembl_gene_id')  # merge conversion result and original dataframe
    merged_res = merged_res.drop(columns=['ensembl_gene_id'])  # we do not need original ensembl gene ID anymore
    merged_res = merged_res.set_index('entrez_gene_id')  # set converted gene ID as a row index

    # handle multiple entries (in case where one ensembl ID mapped to multiple entrez ID)
    merged_res = merged_res.groupby(merged_res.index).mean().round()  # set averaged and rounded-down value for a duplicated rows
    merged_res.sort_index(axis=1, inplace=True)  # sort column names (sample names) alphabetically

    return merged_res
