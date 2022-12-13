import os
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