import os
from typing import Union
import numpy as np
import pandas as pd
from statsmodels.robust import mad




class STARLogAnalyzer:
    def __init__(self, file_path: str) -> None:
        self.file_path = file_path
        self.full_log = self.MultiLogParser()


    def LogParser(self) -> pd.DataFrame:
        """
            Parse Single Log file
            Helper function for MultiLogParser()
        """
        row_names = self.file_path.split('/')[-1][:-14]
    
        with open(self.file_path, 'r') as log:
            cols = []
            rows = []
            
            idx = 0
            while True:
                line = log.readline().strip()  # remove heading and tailing whitespaces
                
                if idx < 5:  # skip first 5 lines (useless info)
                    idx += 1
                    continue
                
                if not line:  # empty line below 5th line means end of file
                    break
                
                if '|' not in line:  # line without '|' means simple informative row
                    continue
                
                col, row = line.split('|')
                cols.append(col.strip())  # remove heading and tailing whitespaces
                rows.append(row.strip())  # remove heading and tailing whitespaces
        
        return pd.DataFrame({c:r for c, r in zip(cols, rows)}, index=[row_names])
    

    def MultiLogParser(self) -> pd.DataFrame:
        """
            Parse Multiple Log files
        """
        res_df = pd.DataFrame()  # empty dataframe
    
        fnames = os.listdir(self.file_path)
        
        for fname in fnames:
            fpath = os.path.join(self.file_path, fname)
            res_df = pd.concat([res_df, self.LogParser(fpath)])

        res_df.sort_index(inplace=True)
        
        return res_df




class CountChef:
    def __init__(self, file_path: str) -> None:
        self.file_path = file_path
        self.cnt, self.fpkm, self.tpm, self.gene_length = self.CounTractor()


    def CounTractor(self) -> tuple:
        cnt_df = pd.DataFrame()
        tpm_df = pd.DataFrame()
        fpkm_df = pd.DataFrame()
        
        rsem_results = os.listdir(self.file_path)
        
        first = True
        for rsem in rsem_results:
            temp = pd.read_csv(os.path.join(self.file_path, rsem), sep='\t')
            temp.index = temp['gene_id']
            
            if first:
                first = False
            else:
                if all(temp.index != cnt_df.index):
                    print("Error, Index not matching!!\n")
                    break
            
            name = rsem[:-14]
            
            cnt_df[name] = temp['expected_count']
            fpkm_df[name] = temp['FPKM']
            tpm_df[name] = temp['TPM']
        
        return (cnt_df, fpkm_df, tpm_df, temp['length'])
    

    def CountQC(self) -> None:
        pass

    
    def MADFilter(self, cnt_df: pd.DataFrame) -> pd.DataFrame:  # this function could be applied to any kinds of numeric dataframes (including TPM, FPKM, ...)
        """
            MADFilter : remove outliers with Median Absolute Deviation
                        each MAD value is calculated per sample (column-wise)
                        MAD filter will be applied per sample
                        final output will have genes that every MAD filter applied sample has
        """
        crit_ub = cnt_df.median(axis=0) + 3*mad(cnt_df, axis=0)  # upper bound : per sample calculation (axis=0)
        crit_lb = cnt_df.median(axis=0) - 3*mad(cnt_df, axis=0)  # lower bound : per sample calculation (axis=0)
        
        # (cnt_df > lower_bound) and (cnt_df < upper_bound) : matrix-wise work --> each element in dataframe results in True or False
        # (True or False dataframe).sum(axis=1) : Gene-wise summation (value of gene A will be number of True samples for gene A)
        # (True or False dataframe).sum(axis=1) == cnt_df.shape[1] : only True if every sample is True for a gene-wise bool operation
        return cnt_df[((cnt_df >= crit_lb) & (cnt_df <= crit_ub)).sum(axis=1) == cnt_df.shape[1]]  # retain genes of which every sample has value of 'True'


    def DataCleansing(self,
                      data: pd.DataFrame,  # this function could be applied to any kinds of numeric dataframes (including TPM, FPKM, ...)
                      hard_min_cnt: int=7,  # 7 = 28(sample number) / 4(group number)
                      soft_min_cnt: int=7,
                      sample_threshold: int=3) -> pd.DataFrame:  # minimum number of samples per group is 3 (same as minimum number of biological replicates)
        """
            DataCleansing : guarantee minimum read counts per sample and apply MAD filter
        """
        #  hard_min_cnt : every sample must have at least har_min_cnt for each gene --> remove gene if any of sample has a count less than hard_min_cnt
        data_cleaned = data[(data >= hard_min_cnt).sum(axis=1) == data.shape[1]]  # retain genes that have counts bigger than hard_min_cnt for every sample
        
        # if you want to give more strict cutoff, you may use soft_min_cnt with manually selected minimum number of satisfying samples
        data_cleaned = data_cleaned[(data_cleaned >= soft_min_cnt).sum(axis=1) >= sample_threshold]  # there should be at least #sample_threshold number of samples that satisfy

        data_cleaned = self.MADFilter(data_cleaned)
        
        return data_cleaned


def SaveObject(object_name: Union[str, list], base_path: str) -> None:

    if isinstance(object_name, list):
        for name in object_name:
            SaveObject(name, base_path)
    
    elif isinstance(object_name, str):
        file_path = os.path.join(base_path, object_name + '.parquet')
    
        if not os.path.exists(file_path):
            print(f"No saved '{object_name}' found. Saving now...!")
            globals()[object_name].to_parquet(file_path)



def LoadObject(file_path: str) -> None:
    for file in os.listdir(file_path):
        name = file.split('.')[0]
        
        if name in globals():
            continue
        
        print("Loading", name, "from parquet data")
        globals().update({name : pd.read_parquet(os.path.join(file_path, file))})