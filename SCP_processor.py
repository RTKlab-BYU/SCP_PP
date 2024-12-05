import os
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
import numpy as np
import pandas as pd
import re
from pandas.api.types import is_numeric_dtype
from scipy.stats import ttest_ind_from_stats
import json

class SCP_processor:
    def __init__(self, ignore_peptides=False, ignore_proteins=False, min_peptides=1, data_type = "LF") -> None:
        self.min_unique_peptides = min_peptides
        self.ignore_peptides = ignore_peptides
        self.ignore_proteins = ignore_proteins
        self.data_type = data_type

    def sumIDs(self,IDMatrix):
        """_summarize the ID matrix infor into ID summary_
        
        Args:
            IDMatrix (_type_): _protein or pepetides matrix_
            0 Accession/Annotated Sequence 	run1 	run2 	run3 
            1 P023D12	MS2 	MBR 	NaN 
            2 P1222	NaN 	ID 	NaN 
        ID: means we don't know the ID mode

        Returns:
            _type_: _description_
                                        names  MS2_IDs  MBR_IDs  Total_IDs
    0            10ng_QC_1_channel2 Intensity      NaN      NaN       3650
    1            10ng_QC_2_channel1 Intensity      NaN      NaN       3604
    ....
            """
        if self.data_type == "TMT":
            # removes the columns that don't have ID data
            columns = [col for col in IDMatrix.columns if not any(
                substring in col for substring in [
                    'Accession', 'Annotated Sequence'])]
            runs = list(set(["-".join(x.split("-")[0:2])  for x in columns]))
            # print(runs)
            #put each ID_Modes into a list
            MS2_ID = []
            MBR_ID = []
            total_ID = []
            for eachRun in runs:
                currentMatrix = IDMatrix.loc[:,IDMatrix.columns.str.contains(eachRun+"-")]
                MS2_ID.append(len(currentMatrix[currentMatrix.isin(["MS2"]).any(axis=1)])) #PD differentiates
                MBR_ID.append(len(currentMatrix[currentMatrix.isin(["MBR"]).any(axis=1)])) #PD differentiates
                total_ID_each = len(currentMatrix[currentMatrix.isin(["ID"]).any(axis=1)]) #some don't so we count total directly
                if total_ID_each == 0: #otherwise we sum
                    total_ID_each = len(currentMatrix[currentMatrix.isin(["MS2"]).any(axis=1)]) + len(currentMatrix[currentMatrix.isin(["MBR"]).any(axis=1)])
                total_ID.append(total_ID_each)
        else:
            runs = [col for col in IDMatrix.columns if not any(substring in col for substring in [
                    'Accession', 'Annotated Sequence'])]
            #put each ID_Modes into a list
            returnNames = []
            MS2_ID = []
            MBR_ID = []
            total_ID = []
            for eachColumn in runs:
                MS2_ID.append(len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "MS2"])) #PD differentiates
                MBR_ID.append(len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "MBR"])) #PD differentiates
                total_ID_each = len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "ID"]) #some don't so we count total directly
                if total_ID_each == 0: #otherwise we sum
                    total_ID_each = len(IDMatrix[eachColumn][
                        IDMatrix[eachColumn] == "MS2"]) + len(IDMatrix[
                        eachColumn][IDMatrix[eachColumn] == "MBR"])
                total_ID.append(total_ID_each)


        return pd.DataFrame({'names': runs,
                            'MS2_IDs': MS2_ID,
                            'MBR_IDs': MBR_ID,
                            'Total_IDs': total_ID})

    def generate_column_from_name_mapping(self,columns, partial_column_name_mapping):
        #input is column names, and a dictionary with what you want each column (key) to be renamed to (value)
        column_name_mapping = {}
        for col in columns:
            for key, value in partial_column_name_mapping.items():
                if str(key) in str(col): #in the case of PD, we are looking for a pattern within the column name
                    column_name_mapping[col] = value
                    break
        return column_name_mapping

    def generate_column_to_name_mapping(self, columns, partial_column_name_mapping):
        """Generates a dictionary that maps column names to a new name
        
        Parameters:
        columns (lst): List of current column names
        partial_column_name_mapping (dict): Map of each column (key) to be renamed to (value)

        Returns
        column_name_mapping (dict): Map of each column (key) and its new name (value)
        """
        column_name_mapping = {}
        for col in columns:
            for key, value in partial_column_name_mapping.items():
                if key in col:  #after we get away from PD's weirdness, then we want exact matches,
                                #so we don't get for example 1-11 when we look for 1-1 for file Identifiers or filenames
                    column_name_mapping[col] = value
                    break
        return column_name_mapping


    def combine_IDs(self,all_matrix, MS2_matrix):
        # make IDs into MBR
        if "Annotated Sequence" in all_matrix.columns:
            name = "Annotated Sequence"
        elif "Accession" in all_matrix.columns:  
            name = "Accession"
        id_cols = all_matrix.columns.tolist()
        id_cols.remove(name)
        
        
        # for eachColumn in id_cols:
        #     if len(all_matrix[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2"))]) > 0:
        #         all_matrix.loc[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2")), [eachColumn]] = MS2_matrix[[eachColumn]]

        all_keys = pd.merge(all_matrix[name], MS2_matrix[name],how="outer")
        all_matrix = pd.merge(all_matrix, all_keys, how="right").replace("ID","MBR")
        MS2_matrix = pd.merge(MS2_matrix, all_keys, how="right")

        for eachColumn in id_cols:
            if len(all_matrix[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2"))]) > 0:
                all_matrix.loc[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2")), [eachColumn]] = "MS2"

        return all_matrix 


    def combine_diann_IDs(self, all_matrix, MS2_matrix):
        # make IDs into MBR
        id_cols = all_matrix.columns.tolist()

        if "Annotated Sequence" in all_matrix.columns:
            name = "Annotated Sequence"
            id_cols.remove("Precursor.Charge")
            all_keys = pd.merge(all_matrix[[name, "Precursor.Charge"]], MS2_matrix[[name, "Precursor.Charge"]], how="outer")
            all_matrix = pd.merge(all_matrix, all_keys, how="right").replace("ID","MBR")
            MS2_matrix = pd.merge(MS2_matrix, all_keys, how="right")

        elif "Accession" in all_matrix.columns:  
            name = "Accession"
            id_cols.remove("Protein.Names")
            all_keys = pd.merge(all_matrix[[name, "Protein.Names"]], MS2_matrix[[name, "Protein.Names"]], how="outer")
            all_matrix = pd.merge(all_matrix, all_keys, how="right").replace("ID","MBR")
            MS2_matrix = pd.merge(MS2_matrix, all_keys, how="right")
    
        id_cols.remove(name)
        for eachColumn in id_cols:
            if len(all_matrix[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2"))]) > 0:
                all_matrix.loc[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2")), [eachColumn]] = "MS2"

        return all_matrix



    def read_file(self,queue_id=None, queue_info=None, processor_info=None,
                input1=None, input2=None, input3=None, input4=None, input5=None,
                process_app=None, file_id=1, USE_MaxLFQ=False, this_organism=None):
        """_Read data from data manager API or through local files or read directly
        in the webapp_
        Args:
            queue_id (_int_): _processing queue id_
            queue_info (_dict_): _queue info from the API_
            processor_info (_dict_): _processor info from the API_        
            input1 (_str_): _input file 1_ 
            input2 (_str_): _input file 2_
            input3 (_str_): _input file 3_
            input4 (_str_): _input file 4_
            input5 (_str_): _input file 5_
            process_app (_str_): _process app name_
        Returns:
            _dict_: _dictionary containing data all data        
        """

        """ Input files are as followws
        App       Input file
                1                         2                       3         4         5
        PD        _Proteins                 _PeptideGroups                              _InputFiles
        Fragpipe  combined_protein          combined_peptide
        DIANN     diann-output.pg_matrix    diann-output.pr_matrix  protein   peptide   filelist_diann.txt
        """

        if "FragPipe_LF" in process_app:     # fragpipe results
            # read data
            if self.ignore_peptides == False:
                peptide_table = pd.read_table(input2,low_memory=False)
                peptide_table= peptide_table[~peptide_table[
                'Protein'].str.contains(
                "contam_sp", na=False)]

                if USE_MaxLFQ:
                    column_tail = " MaxLFQ Intensity"
                    intensity_cols = peptide_table.columns[
                        peptide_table.columns.str.contains(column_tail)].tolist() #These should be the same for MS2
                else:
                    column_tail = " Intensity" #this includes MaxLFQ Intensity
                    cols = peptide_table.columns
                    intensity_cols =  peptide_table.columns[( \
                        cols.str.contains(" Intensity")) & \
                        (~cols.str.contains(" MaxLFQ Intensity"))]
                cols = peptide_table.columns
                pep_ID_cols = peptide_table.columns[(cols.str.contains(" Spectral Count"))] #so we remove MaxLFQ Intensity here
                #get the column names of peptide_table that contain "Intensity" but not " MaxLFQ"
                ## Peptide abundance table
                peptide_table.rename(columns={'Peptide Sequence': 'Annotated Sequence'}, inplace=True)
                all_abundance_cols = intensity_cols.append(pd.Index(['Annotated Sequence']))
                pep_abundance = peptide_table.loc[:, all_abundance_cols]

                # Proteins ID table
                all_ID_cols = pep_ID_cols.append(pd.Index(['Annotated Sequence']))
                pep_ID_MS2 = peptide_table.loc[:, all_ID_cols]

                pep_abundance = pep_abundance.rename(columns={
                col: col.replace(column_tail, "") for col in
                pep_abundance.columns if column_tail in col})
            
           

                pep_ID_MS2 = pep_ID_MS2.rename(columns={
                    col: col.replace(" Spectral Count", "") for col in
                    pep_ID_MS2.columns if " Spectral Count" in col})
                # remove "Spectral Count", " MaxLFQ Intensity" or " Intensity" from names
                run_name_list = [name.replace(column_tail, "")
                                    for name in intensity_cols]
                
                run_name_list = pd.DataFrame({"Run Names": run_name_list})
                run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

                for item in [pep_abundance,pep_ID_MS2]:
                # Generate a new column name mapping using the function
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
                pep_ID = pep_abundance.copy()
                cols = [col for col in pep_ID.columns if col != 'Annotated Sequence	']
                for col in cols:
                    if pep_ID[col].dtype != 'object': # Check if not a string column
                        pep_ID[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        pep_ID[col] = pep_ID[col].astype(str).str.replace(r"\d+\.\d+", "ID", regex=True)
                #Rename protein numbers
                cols = [col for col in pep_ID_MS2.columns if col != 'Annotated Sequence	']
                for col in cols:
                    if pep_ID_MS2[col].dtype != 'object': # Check if not a string column
                        pep_ID_MS2[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace(r"\d+\.\d+", "MS2", regex=True)

                # print(run_name_list)
                pep_ID = self.combine_IDs(pep_ID, pep_ID_MS2)
                
                pep_other_info = peptide_table.loc[
                    :, ~peptide_table.columns.str.contains('Intensity')]
                pep_other_info["Source_File"] = input2

            if self.ignore_proteins == False:
                protein_table = pd.read_table(input1,low_memory=False)

                # filter Contaminant
                
                protein_table= protein_table[~protein_table['Protein'].str.contains(
                    "contam_sp", na=False)].query(
                    "`Combined Total Peptides` >= @self.min_unique_peptides")
                # get experiment names from columns names containning "Intensity"
                # or " MaxLFQ Intensity" if MaxLFQ is used

                if USE_MaxLFQ:
                    column_tail = " MaxLFQ Intensity"
                    intensity_cols = protein_table.columns[
                        protein_table.columns.str.contains(column_tail)].tolist() #These should be the same for MS2
                else:
                    column_tail = " Intensity" #this includes MaxLFQ Intensity
                    cols = protein_table.columns
                    intensity_cols =  protein_table.columns[( \
                        cols.str.contains(" Intensity")) & \
                        (~cols.str.contains(" MaxLFQ Intensity"))]
                prot_ID_cols = protein_table.columns[(cols.str.contains(" Unique Spectral Count")) & ~(cols == "Combined Unique Spectral Count")] #so we remove MaxLFQ Intensity here
            
            
                # ALL
                ## Proteins abundance table
                protein_table.rename(columns={'Protein ID': 'Accession'},inplace=True)
                all_abundance_cols = intensity_cols.append(pd.Index(['Accession']))
                prot_abundance = protein_table.loc[:, all_abundance_cols]
            

                # Proteins ID table
                all_ID_cols = prot_ID_cols.append(pd.Index(['Accession']))
                prot_ID_MS2 = protein_table.loc[:, all_ID_cols]
                prot_abundance = prot_abundance.rename(columns={
                col: col.replace(column_tail, "") for col in
                prot_abundance.columns if column_tail in col})
                prot_ID_MS2 = prot_ID_MS2.rename(columns={
                col: col.replace(" Unique Spectral Count", "") for col in
                prot_ID_MS2.columns if " Unique Spectral Count" in col})

                # remove "Spectral Count", " MaxLFQ Intensity" or " Intensity" from names
                run_name_list = [name.replace(column_tail, "")
                                    for name in intensity_cols]
                
                run_name_list = pd.DataFrame({"Run Names": run_name_list})
                run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

                for item in [prot_abundance,prot_ID_MS2]:
                # Generate a new column name mapping using the function
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
            

                # get ID matrix tables
                prot_ID = prot_abundance.copy()
                cols = [col for col in prot_ID.columns if col != 'Accession']
                for col in cols:
                    if prot_ID[col].dtype != 'object': # Check if not a string column
                        prot_ID[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        prot_ID[col] = prot_ID[col].astype(str).str.replace(r"\d+\.\d+", "ID", regex=True)
                            
            
                cols = [col for col in prot_ID_MS2.columns if col != 'Accession']
                for col in cols:
                    if prot_ID_MS2[col].dtype != 'object': # Check if not a string column
                        prot_ID_MS2[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace(r"\d+\.\d+", "MS2", regex=True)
                
                prot_ID = self.combine_IDs(prot_ID, prot_ID_MS2)

                prot_other_info = protein_table.loc[
                    :, ~protein_table.columns.str.contains('Intensity')]
                prot_other_info["Source_File"] = input1
            

        elif "PD_LF" in process_app:

            meta_table = pd.read_table(input5,low_memory=False)
            #filter rows in meta table on File ID column if it is NaN
            meta_table = meta_table[meta_table['File ID'].notna()]
            if self.ignore_peptides == False:
                pep_other_info["Source_File"] = input2
                peptide_table = pd.read_table(input2,low_memory=False)  
                peptide_table= peptide_table[(peptide_table[
                'Contaminant'] == False) & (peptide_table["Confidence"]== "High")]

                pep_abundance = peptide_table.filter(regex='Abundance:|Annotated Sequence')
                pep_ID = peptide_table.filter(regex='Found in Sample:|Annotated Sequence')
                pep_other_info = peptide_table.loc[:, ~peptide_table.columns.str.contains('Found in Sample:|Abundance:')]
                new_dict = {"Abundance: " + key + ":": value for key, value in file_path_name_dict.items()}
                for item in [pep_abundance]:
                    # Generate a new column name mapping using the function
                    column_name_mapping = self.generate_column_from_name_mapping(item.columns, new_dict)
                    #TODO solving  A value is trying to be set on a copy of a slice from a DataFrame
                    
                    item.rename(columns = column_name_mapping, inplace = True)
                    #use generate_column_to_name_mapping function because we don't want partial matches as in QC_HeLa.raw and QC_HeLa_20230727235101.raw
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
                

                new_dict = {"Found in Sample: " + key + ":": value for key, value in file_path_name_dict.items()}
                for item in [pep_ID]:
                    # Generate a new column name mapping using the function
                    column_name_mapping = self.generate_column_from_name_mapping(item.columns, new_dict)

                    item.rename(columns = column_name_mapping, inplace = True)
                    #use generate_column_to_name_mapping function because we don't want partial matches
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))

                    item.rename(columns = fileid_mapping,inplace=True)

                # replace "High" to MS2 "Peak Found" to MBR, the rest become np.NaN
                replacements = {'High': 'MS2', 'Peak Found': 'MBR', "Medium": np.NaN, "Low": np.NaN, "Not Found": np.NaN}
                for column in run_name_list["Run Identifier"]:
                    if column in pep_ID.columns:
                        pep_ID[column] = pep_ID[column].replace(to_replace=replacements)
            
                peptide_ID_summary = self.sumIDs(pep_ID)
            if self.ignore_proteins == False:
                protein_table = pd.read_table(input1,low_memory=False)
                
                # filter Contaminant
                protein_table= protein_table[(protein_table[
                    "Protein FDR Confidence: Combined"] == "High") &
                                ((protein_table["Master"] == "IsMasterProtein") | 
                                (protein_table["Master"] == "Master")) & 
                                (protein_table["Contaminant"] == False)]

                protein_table.rename(
                    columns={'# Peptides': 'number of peptides'}, inplace=True)
                protein_table=protein_table.query(
                    "`number of peptides` >= @self.min_unique_peptides")
                prot_abundance = protein_table.filter(regex='Abundance:|Accession')
                prot_ID = protein_table.filter(regex='Found in Sample:|Accession')
                prot_other_info = protein_table.loc[:, ~protein_table.columns.str.contains('Found in Sample:|Abundance:')]
                #format the read in table into three different tables: abundance, id and other_info
            
                prot_other_info["Source_File"] = input1


            #change column names to file/run names to our fileID

                new_dict = {"Abundance: " + key + ":": value for key, value in file_path_name_dict.items()}
                for item in [prot_abundance]:
                    # Generate a new column name mapping using the function
                    column_name_mapping = self.generate_column_from_name_mapping(item.columns, new_dict)
                    #TODO solving  A value is trying to be set on a copy of a slice from a DataFrame
                    
                    item.rename(columns = column_name_mapping, inplace = True)
                    #use generate_column_to_name_mapping function because we don't want partial matches as in QC_HeLa.raw and QC_HeLa_20230727235101.raw
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
                

                new_dict = {"Found in Sample: " + key + ":": value for key, value in file_path_name_dict.items()}
                for item in [prot_ID]:
                    # Generate a new column name mapping using the function
                    column_name_mapping = self.generate_column_from_name_mapping(item.columns, new_dict)

                    item.rename(columns = column_name_mapping, inplace = True)
                    #use generate_column_to_name_mapping function because we don't want partial matches
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))

                    item.rename(columns = fileid_mapping,inplace=True)

                # replace "High" to MS2 "Peak Found" to MBR, the rest become np.NaN
                replacements = {'High': 'MS2', 'Peak Found': 'MBR', "Medium": np.NaN, "Low": np.NaN, "Not Found": np.NaN}
                for column in run_name_list["Run Identifier"]:
                    
                    if column in prot_ID.columns:
                        prot_ID[column] = prot_ID[column].replace(to_replace=replacements)
            
                # get ID summary by parsing ID Matrix
                protein_ID_summary = self.sumIDs(prot_ID)
            

            # Replace single backslashes with forward slashes in the 'file_paths' column
            meta_table['File Name'] = meta_table['File Name'].str.replace('\\', '/', regex=False)
            # Apply a lambda function to extract file names without extensions
            meta_table['file_names'] = meta_table['File Name'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
            file_path_name_dict = dict(zip(meta_table['File ID'], meta_table['file_names']))
            run_name_list = pd.DataFrame({"Run Names": file_path_name_dict.values()})
            run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))
            
        elif "DIANN_LF" in process_app:
            ## Peptides
            if self.ignore_peptides is False:

                # read in DIANN output files
                peptide_table = pd.read_table(input2,low_memory=False)
                if this_organism != None:
                    peptide_table = peptide_table.loc[peptide_table["Protein.Names"].str.endswith("_"+this_organism)]
                # We currently do not have a MS2 file
                peptide_table_MS2 = pd.read_table(input4,low_memory=False)
                if this_organism != None:
                    peptide_table_MS2 = peptide_table_MS2.loc[peptide_table_MS2["Protein.Names"].str.endswith("_"+this_organism)]

                # Creates a new data frame
                pep_other_info = pd.DataFrame({"Mapped Proteins": peptide_table["Protein.Group"], "Annotated Sequence": peptide_table["Modified.Sequence"]})
                pep_other_info["Source_File"] = "None"

                # Filters out proteins with "contam" or that are null
                peptide_table= peptide_table[~peptide_table['Protein.Group'].str.contains("contam", na=False)]
                peptide_table_MS2= peptide_table_MS2[~peptide_table_MS2['Protein.Group'].str.contains("contam", na=False)]
                pep_other_info= pep_other_info[~pep_other_info['Mapped Proteins'].str.contains("contam", na=False)]

                # I changed the name in the creation of pep_other_info
                # pep_other_info.rename(columns={'Modified.Sequence': 'Annotated Sequence'}, inplace=True)
                
                # # Replace backslashes with forward slashes if data comes from Windows
                # meta_table['File Name'] = meta_table['File Name'].str.replace('\\', '/', regex=False)
                # # Apply a lambda function to extract file names without extensions
                # meta_table['File Name'] = meta_table['File Name'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
                # run_name_list = meta_table['File Name'].tolist()
                # run_name_list = pd.DataFrame({"Run Names": run_name_list})
                # run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))
               
                # Separate the columns of just "Modified.Sequence" and file names
                peptide_path_cols = peptide_table_MS2.filter(regex='\\\\|Modified.Sequence|Precursor.Charge').columns

                # Create data frames that are just the columns with "Modified.Sequence" or a filename
                pep_abundance = peptide_table.loc[:, peptide_path_cols]
                pep_abundance_MS2 = peptide_table_MS2.loc[:, peptide_path_cols]

                # Rename the column "Modified.Sequence" to "Annotated Sequence"
                if this_organism != None:
                    pep_abundance.columns = [os.path.splitext(os.path.basename(x))[0]+"_"+this_organism if x in file_path_cols else x for x in pep_abundance.columns]
                pep_abundance = pep_abundance.rename(columns={'Modified.Sequence': 'Annotated Sequence'})
                if this_organism != None:
                    pep_abundance_MS2.columns = [os.path.splitext(os.path.basename(x))[0]+"_"+this_organism if x in file_path_cols else x for x in pep_abundance_MS2.columns]
                pep_abundance_MS2 = pep_abundance_MS2.rename(columns={'Modified.Sequence': 'Annotated Sequence'})

                # Create a list of just the file names
                file_path_cols = peptide_table.filter(regex='\\\\').columns

                # Creates a data frame of just the file names (no .raw)
                run_name_list = pd.DataFrame(data={"Run Names": [os.path.splitext(x)[0] for x in file_path_cols]})

                # Creates a new column in the data frame that is used to identify the run by the index of the data frame
                run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

                # Rename the file columns of pep_abundance to be just the run names (no .raw) 
                for item in [pep_abundance,pep_abundance_MS2]:
                    # Yes it is reduntant to make partial_column_name_mapping, but it is to use the function generate_col...
                    filename_to_identifier_map = dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"]))
                    # filename_to_filename_map = {name: name for name in run_name_list["Run Names"]}
                    run_name_mapping = self.generate_column_to_name_mapping(item.columns, filename_to_identifier_map)
                    item.rename(columns = run_name_mapping, inplace=True)

                # Analyze peptide ID
                # Create a copy of pep_abundance to manipulate
                pep_ID = pep_abundance.copy()
                pep_abundance.drop(["Precursor.Charge"], axis=1, inplace=True)

                # Create a list of columns, but ignore the "Annotated Sequence" column
                cols = [col for col in pep_ID.columns if col != 'Annotated Sequence' and col != "Precursor.Charge"]

                # For all of the columns
                for col in cols:
                    # If the column is not an "object" type (if it is not a string)
                    if pep_ID[col].dtype != 'object':
                        # Replace all 0's with NaN
                        pep_ID[col] = pep_ID[col].replace(0, np.nan)
                        # Replace all numerical values (floats, ints, and exponentials) to ID
                        pep_ID[col] = pep_ID[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "ID", regex=True)

                # Repeat the process above with MS2
                pep_ID_MS2 = pep_abundance_MS2.copy()
                pep_abundance_MS2.drop(["Precursor.Charge"], axis=1, inplace=True)
                cols = [col for col in pep_ID_MS2.columns if col != 'Annotated Sequence' and col != "Precursor.Charge"]
                for col in cols:
                    if pep_ID_MS2[col].dtype != 'object': 
                        pep_ID_MS2[col] = pep_ID_MS2[col].replace(0, np.nan)
                        pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "MS2", regex=True)
                
                # Returns a combination of the MS2 and normal data.
                # There is a "MS2" when a run has data for a peptide in the MS2 
                # There is an "MBR" when a run has data for a peptide in just the normal data
                pep_ID = self.combine_diann_IDs(pep_ID, pep_ID_MS2)
                pep_ID.drop(["Precursor.Charge"], axis=1, inplace=True)


                
            ## Proteins
            if self.ignore_proteins is False:

                # Read in DIANN output files (normal)
                protein_table = pd.read_table(input1,low_memory=False)
                if this_organism != None:
                    protein_table = protein_table.loc[protein_table["Protein.Names"].str.endswith("_"+this_organism)]

                # Read in DIANN output files (normal)
                protein_table_MS2 = pd.read_table(input3,low_memory=False)
                if this_organism != None:
                    protein_table_MS2 = protein_table_MS2.loc[protein_table_MS2["Protein.Names"].str.endswith("_"+this_organism)]
                    
                # Create a new data frame with just the Genes and Protein.Group columns
                prot_other_info = pd.DataFrame({"Accession": protein_table["Genes"], "Protein.Group": protein_table["Protein.Group"]})
                prot_other_info["Source_File"] = "None"

                # meta_table = pd.read_csv(input5, sep=' ', header=None, names=["File Name"])

                # filter Contaminant
                protein_table= protein_table[~protein_table['Protein.Group'].str.contains("contam", na=False)]
                protein_table_MS2= protein_table_MS2[~protein_table_MS2['Protein.Group'].str.contains("contam", na=False)]
                prot_other_info= prot_other_info[~prot_other_info['Accession'].str.contains("contam", na=False)]
            
                # Alex took out, changed name in creation of prot_other_info
                # prot_other_info.rename(columns={'Protein': 'Accession'}, inplace=True)

                # Get the file names from the meta table
                protein_path_cols = protein_table_MS2.filter(regex='\\\\|Genes|Protein.Names').columns

                # Create data frames that are just the columns with "Genes" or a filename
                prot_abundance = protein_table.loc[:, protein_path_cols]
                prot_abundance_MS2 = protein_table_MS2.loc[:, protein_path_cols]

                # Create a list of just the file names
                file_path_cols = protein_table.filter(regex='\\\\').columns

                # Changes the columns of prot_abundance to the path/filename (no .raw)
                prot_abundance.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance.columns]
                prot_abundance_MS2.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance_MS2.columns]

                # Rename Genes to Accession, and removes ";----" from the column
                prot_abundance = prot_abundance.rename(columns={'Genes': 'Accession'})
                prot_abundance["Accession"] =  prot_abundance["Accession"].str.replace(";.*","",regex = True)
                prot_abundance_MS2 = prot_abundance_MS2.rename(columns={'Genes': 'Accession'})   
                prot_abundance_MS2["Accession"] =  prot_abundance_MS2["Accession"].str.replace(";.*","",regex = True)

                # Create a new data frame with just the file names and the run identifier "0-n"
                run_name_list = pd.DataFrame(data={"Run Names": [os.path.splitext(os.path.basename(x))[0] for x in file_path_cols]})
                run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

                for item in [prot_abundance, prot_abundance_MS2]:
                    # Generate a new column name mapping using the function
                    filename_to_identifier_map = dict(zip(run_name_list["Run Names"], run_name_list["Run Identifier"]))
                    # filename_to_filename_map = {name: name for name in run_name_list["Run Names"]}
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, filename_to_identifier_map)
                    item.rename(columns=fileid_mapping, inplace=True)


            # Analyze Protein ID
                # Create a copy of prot_abundance to manipulate
                prot_ID = prot_abundance.copy()
                prot_abundance.drop(["Protein.Names"], axis=1, inplace=True)

                # Create a list of columns, but ignore the "Accession" column
                cols = [col for col in prot_ID.columns if col != 'Accession' and col != "Protein.Names"]

                # For all of the columns
                for col in cols:
                    # If the column is not an "object" type (if it is not a string)
                    if prot_ID[col].dtype != 'object':
                        # Replace all 0's with NaN
                        prot_ID[col] = prot_ID[col].replace(0, np.nan)
                        # Replace all numerical values (floats, ints, and exponentials) to ID
                        prot_ID[col] = prot_ID[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "ID", regex=True)

                # Repeat the the process above with MS2
                prot_ID_MS2 = prot_abundance_MS2.copy()
                prot_abundance_MS2.drop(["Protein.Names"], axis=1, inplace=True)
                cols = [col for col in prot_ID_MS2.columns if col != 'Accession' and col != "Protein.Names"]
                for col in cols:
                    if prot_ID_MS2[col].dtype != 'object':
                        prot_ID_MS2[col] = prot_ID_MS2[col].replace(0, np.nan)
                        prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "MS2", regex=True)


                # Returns a combination of the MS2 and normal data.
                # There is a "MS2" when a run has data for a peptide in the MS2 
                # There is an "MBR" when a run has data for a peptide in just the normal data
                prot_ID = self.combine_diann_IDs(prot_ID, prot_ID_MS2)   
                prot_ID.drop(["Protein.Names"], axis=1, inplace=True)

            
        elif "FragPipe_TMT" in process_app:     # fragpipe results
            # read data
            if self.ignore_peptides == False:
                peptide_table = pd.read_table(input2,low_memory=False)
                #remove 
                #protein remove = 

                # pep_info_columns = ["Index", "Gene", "ProteinID", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity"]
                pep_info_columns = ["Index", "Gene", "ProteinID", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity"]
                 ## Peptide abundance table
                peptide_table.rename(columns={'Peptide': 'Annotated Sequence'}, inplace=True)
                pep_abundance = peptide_table.drop(pep_info_columns,axis=1)
                pep_info_columns.append("Annotated Sequence")
                pep_other_info = peptide_table.loc[:,pep_info_columns]
                all_ID_cols = pep_abundance.columns
                pep_ID_MS2 = peptide_table.loc[:, all_ID_cols]
                run_list = [name.split("#")[0] for name in all_ID_cols.drop("Annotated Sequence")]
                channel_list = [name.split("#")[1] for name in all_ID_cols.drop("Annotated Sequence")]
                
                run_name_list = pd.DataFrame({"Run Names": run_list})
                unique_run_list = list(set(run_list))
                run_ids = []
                channel_ids = []
                for run in range(len(unique_run_list)):
                    all_channels = run_name_list.loc[run_name_list['Run Names'] == unique_run_list[run]]
                    run_ids = run_ids+ ["-".join([str(file_id),str(run)]) for x in range(len(all_channels))]
                    channel_ids = channel_ids + ["-".join([str(file_id),str(run),str(i)]) for i in range(len(all_channels))]
                # print(pep_ID_MS2.columns)
                run_name_list["Run Identifier"] = run_ids
                run_name_list["Channel Identifier"] = channel_ids


                for item in [pep_abundance,pep_ID_MS2]:
                    # Generate a new column name mapping using the function
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(all_ID_cols.drop("Annotated Sequence"),run_name_list["Channel Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
                

                # get ID matrix tables
                pep_ID = pep_abundance.copy()
                cols = [col for col in pep_ID.columns if col != 'Annotated Sequence	']
                for col in cols:
                    if pep_ID[col].dtype != 'object': # Check if not a string column
                        pep_ID[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        pep_ID[col] = pep_ID[col].astype(str).str.replace(r"\d+\.\d+", "ID", regex=True)
                cols = [col for col in pep_ID_MS2.columns if col != 'Annotated Sequence	']
                for col in cols:
                    if pep_ID_MS2[col].dtype != 'object': # Check if not a string column
                        pep_ID_MS2[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace(r"\d+\.\d+", "MS2", regex=True)
                pep_ID = self.combine_IDs(pep_ID, pep_ID_MS2)
                pep_other_info["Source_File"] = input2  

            if self.ignore_proteins == False:
                protein_table = pd.read_table(input1,low_memory=False)
                
            
                prot_info_columns = ["Index",  "NumberPSM", "MaxPepProb", "ReferenceIntensity"]

                # ALL
                ## Proteins abundance table
                protein_table.rename(columns={'Gene': 'Accession'},inplace=True)
                prot_abundance = protein_table.drop(prot_info_columns,axis=1)
                prot_info_columns.append("Accession")
                prot_other_info = protein_table.loc[:,prot_info_columns]
                
            
            # Proteins ID table
                all_ID_cols = prot_abundance.columns
                prot_ID_MS2 = protein_table.loc[:, all_ID_cols]
                
            # Proteins ID table
            
            # remove "Spectral Count", " MaxLFQ Intensity" or " Intensity" from names
                run_list = [name.split("#")[0] for name in all_ID_cols.drop("Accession")]
                channel_list = [name.split("#")[1] for name in all_ID_cols.drop("Accession")]
                
                run_name_list = pd.DataFrame({"Run Names": run_list})
                unique_run_list = list(set(run_list))
                run_ids = []
                channel_ids = []
                for run in range(len(unique_run_list)):
                    all_channels = run_name_list.loc[run_name_list['Run Names'] == unique_run_list[run]]
                    run_ids = run_ids+ ["-".join([str(file_id),str(run)]) for x in range(len(all_channels))]
                    channel_ids = channel_ids + ["-".join([str(file_id),str(run),str(i)]) for i in range(len(all_channels))]
                # print(pep_ID_MS2.columns)
                run_name_list["Run Identifier"] = run_ids
                run_name_list["Channel Identifier"] = channel_ids


                for item in [prot_abundance,prot_ID_MS2]:
                    # Generate a new column name mapping using the function
                    fileid_mapping = self.generate_column_to_name_mapping(item.columns, dict(zip(all_ID_cols.drop("Accession"),run_name_list["Channel Identifier"])))
                    item.rename(columns = fileid_mapping,inplace=True)
                

                # get ID matrix tables
                prot_ID = prot_abundance.copy()
                cols = [col for col in prot_ID.columns if col != 'Accession']
                for col in cols:
                    if prot_ID[col].dtype != 'object': # Check if not a string column
                        prot_ID[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        prot_ID[col] = prot_ID[col].astype(str).str.replace(r"\d+\.\d+", "ID", regex=True)
                pep_ID = pep_abundance.copy()
            #Rename protein numbers
            
                cols = [col for col in prot_ID_MS2.columns if col != 'Accession']
                for col in cols:
                    if prot_ID_MS2[col].dtype != 'object': # Check if not a string column
                        prot_ID_MS2[col].replace(0, np.nan, inplace=True)
                        # Replace all numerical values to ID
                        prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace(r"\d+\.\d+", "MS2", regex=True)
            
            # print(run_name_list)
                prot_ID = self.combine_IDs(prot_ID, prot_ID_MS2)    

                prot_other_info["Source_File"] = input1
               
        # elif "MQ_TMT" in process_app:
        #     protein_table = pd.read_table(input1,low_memory=False)
        #     peptide_table = pd.read_table(input2,low_memory=False)
            
        #     # filter Contaminant
        #     protein_table= protein_table[
        #                     (protein_table["Potential contaminant"] != "+")]
        #     protein_table.rename(
        #         columns={'Unique peptides': 'number of peptides'}, inplace=True)
        #     protein_table = protein_table[(protein_table["number of peptides"] > self.min_unique_peptides)]
        #     peptide_table= peptide_table[(peptide_table[
        #         'Potential contaminant'] != "+") ]

        #     meta_table = pd.read_table(input5,low_memory=False)
        #     #filter out last row in meta table because it is Total
        #     meta_table = meta_table[:-1]

           
        #     run_name_list = pd.DataFrame({"Run Names": meta_table["Raw file"]})
        #     run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))
            
        #     #format the read in table into three different tables: abundance, id and other_info
        #     prot_abundance = protein_table.filter(regex='Reporter intensity |Gene names')
        #     prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity corrected ")]
        #     prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity count ")]
        #     prot_ID = protein_table.filter(regex='Identification type |Gene names')
        #     prot_other_info = protein_table.loc[:, ~protein_table.columns.str.contains('Identification type |Reporter intensity ')]
            

        #     pep_abundance = peptide_table.filter(regex='Reporter intensity |Sequence')
        #     pep_abundance = pep_abundance.loc[:,~pep_abundance.columns.str.contains("Reporter intensity corrected ")]
        #     pep_abundance = pep_abundance.loc[:,~pep_abundance.columns.str.contains("Reporter intensity count ")]
        #     pep_ID = peptide_table.filter(regex='Identification type |Sequence')
        #     pep_other_info = peptide_table.loc[:, ~peptide_table.columns.str.contains('Identification type |Reporter intensity ')]

        #     prot_other_info["Source_File"] = input1
        #     pep_other_info["Source_File"] = input2
            
        #     #change column names to file/run names to our fileID
        #     temp_table = meta_table.copy()
        #     temp_table["start"] = "Identification type "
        #     if is_numeric_dtype(meta_table["Experiment"]):
        #         temp_table["Experiment"] = temp_table["Experiment"].astype(str).str.replace('\.[0-9]+', '',regex=True)
        #     else:
        #         print("*********")
        #     temp_table["pattern"] = temp_table["start"]+temp_table["Experiment"]
        #     new_dict =  dict(zip(temp_table["pattern"] ,run_name_list["Run Identifier"]))
           
        #     new_dict["Gene names"] = "Accession"
        #     new_dict["Sequence"] = "Annotated Sequence"
        #     for item in [pep_ID,prot_ID,prot_other_info,pep_other_info]:
        #         # Generate a new column name mapping using the function
        #         fileid_mapping = self.generate_column_to_name_mapping(item.columns,new_dict)
        #         item.rename(columns = fileid_mapping,inplace=True)

        #     #change column names to file/run names to our fileID
        #     old_dict = dict(zip(temp_table["Experiment"] ,run_name_list["Run Identifier"]))
        #     old_dict["Gene names"] = "Accession"
        #     old_dict["Sequence"] = "Annotated Sequence"
        #     max_channel = 0 
        #     for item in [prot_abundance,pep_abundance]:
        #         # Generate a new column name mapping using the function                
        #         i =0
        #         for eachcol in item.columns.tolist():
        #             if eachcol not in old_dict.keys():
        #                 new_dict[eachcol]=old_dict[" ".join(eachcol.split(" ")[3:])] +"-"+eachcol.split(" ")[2]
        #                 max_channel = max(max_channel,int(eachcol.split(" ")[2]))
        #             i = i + 1                
        #             item.loc[item[eachcol]==0,eachcol]=np.NaN
        #         fileid_mapping = self.generate_column_to_name_mapping(item.columns,new_dict)

        #         item.rename(columns = fileid_mapping,inplace=True)
        #         # print(item.columns)
        #     #use generate_column_to_name_mapping function because we don't want partial matches
        #     # replace "High" to MS2 "Peak Found" to MBR, the rest become np.NaN
        #     replacements = {'By MS/MS': 'MS2', 'By matching': 'MBR',  "None": np.NaN}

        #     new_run_name_list = pd.DataFrame({"Run Names":[],"Run Identifier":[]})
        #     for index, each_row in run_name_list.iterrows():
        #         old_identifier = each_row["Run Identifier"]
        #         for channel in range(1,max_channel+1):
        #             each_row["Channel Identifier"] =  old_identifier + "-" + str(channel)
        #             new_run_name_list = pd.concat([new_run_name_list,each_row.to_frame().T])
                
        #     run_name_list = new_run_name_list

        #     for column in run_name_list["Run Identifier"].drop_duplicates():
        #         if column in pep_ID.columns:
        #             for channel in run_name_list.loc[run_name_list["Run Identifier"]==column,"Channel Identifier"].tolist():
        #                 pep_ID[channel] = pep_ID[column]
        #                 pep_ID[channel] = pep_ID[channel].replace(to_replace=replacements)
        #                 pep_abundance.loc[pd.isna(pep_ID[channel]),channel] = np.NaN
        #             pep_ID = pep_ID.drop(columns=column)
        #         if column in prot_ID.columns:
        #             for channel in run_name_list.loc[run_name_list["Run Identifier"]==column,"Channel Identifier"].tolist():
        #                 # print(prot_abundance.loc[prot_abundance[channel]==np.NaN,channel])
        #                 # print(prot_abundance[channel])
        #                 prot_ID[channel] = prot_ID[column]
        #                 prot_ID[channel] = prot_ID[channel].replace(to_replace=replacements)
        #                 prot_abundance.loc[pd.isna(prot_ID[channel]),channel] = np.NaN
        #                 # print(prot_abundance.loc[pd.isna(prot_ID[channel]),channel])
        #             prot_ID = prot_ID.drop(columns=column)
        #     #remove the fake ids (all zeros)
            
        #     prot_abundance = prot_abundance.set_index("Accession")
        #     pep_abundance = pep_abundance.set_index("Annotated Sequence")
        #     for run in run_name_list["Run Identifier"].drop_duplicates():
        #         current_channels = run_name_list.loc[run_name_list["Run Identifier"]==run,"Channel Identifier"].tolist()
        #         # print(current_channels)
        #         # print(prot_ID.columns)
        #         # Identify rows in df2 where all values in the subset are NaN
        #         rows_with_all_nan = prot_abundance[current_channels].isna().all(axis=1)

        #         # Extract the IDs of the rows where all values in the subset are NaN
        #         nan_ids = prot_abundance[rows_with_all_nan].index

        #         # Replace corresponding rows in df1 for the same subset of columns with NaN
        #         prot_ID.loc[prot_ID['Accession'].isin(nan_ids), prot_ID.columns.isin(current_channels)] = np.NaN

        #         rows_with_all_nan = pep_abundance[current_channels].isna().all(axis=1)

        #         # Extract the IDs of the rows where all values in the subset are NaN
        #         nan_ids = pep_abundance[rows_with_all_nan].index

        #         # Replace corresponding rows in df1 for the same subset of columns with NaN
        #         pep_ID.loc[pep_ID['Annotated Sequence'].isin(nan_ids), pep_ID.columns.isin(current_channels)] = np.NaN

        #     pep_abundance = pep_abundance.reset_index()
        #     prot_abundance = prot_abundance.reset_index()

        #     # get ID summary by parsing ID Matrix
        #     protein_ID_summary = self.sumIDs(prot_ID)
        #     peptide_ID_summary = self.sumIDs(pep_ID)
        
        
        # Creates a column in run_name_list for the processing app 
        run_name_list["Processing App"] = process_app

        return_matrix = {'run_metadata': run_name_list}  
        if self.ignore_peptides is False:
            peptide_ID_summary = self.sumIDs(pep_ID)  
            return_matrix['peptide_other_info']=  pep_other_info
            return_matrix['peptide_abundance']=  pep_abundance
            return_matrix['peptide_ID_matrix']=  pep_ID
            return_matrix['peptide_ID_Summary']=  peptide_ID_summary

        if self.ignore_proteins is False:
            protein_ID_summary = self.sumIDs(prot_ID)
            return_matrix['protein_other_info']= prot_other_info
            return_matrix['protein_abundance']=  prot_abundance
            return_matrix['protein_ID_matrix']=  prot_ID
            return_matrix['protein_ID_Summary']=  protein_ID_summary

        return return_matrix

    def read_files(self, queue_ids=None, queue_info=None, processor_info=None, 
                   grouped_input_files=[], organisms=[None]):
        '''
        Creates a list of data objects
        
        Input
        grouped_input_files
        [
            #File 0
        {input1:
        input2:
        input3:   
        input4:
        input5:
        process_app:
        },
            #File 1
        {input1:
        input2:
        input3:   
        input4:
        input5:
        process_app:
        },
        ...
        ]
        '''

        data_objects = []

        i = 0
        for eachGroup in grouped_input_files:
            if queue_ids is not None:
                pass
            else:
                process_app = eachGroup["process_app"]
                input1= eachGroup["input1"]
                input2= eachGroup["input2"]  
                input3= eachGroup["input3"]
                input4= eachGroup["input4"]  
                input5= eachGroup["input5"]

            for each_organism in organisms:
                current_data_object = self.read_file(this_organism=each_organism, input1=input1, input2=input2,
                                                    input3=input3, input4=input4, input5=input5, 
                                                    process_app=process_app, file_id=i)
                data_objects.append(current_data_object)
                i = i + 1

        return data_objects


    def outer_join_data_objects(self,data_objects):
        '''
        Takes in a list of data objects as given by read_files and converts them to a single data object as given by read_files,
        protein info continues to show what was found on each original file, and so forth.
        '''

        first_file = True
        for eachDataObject in data_objects:
            print("***")
            if first_file:
                first_file = False
                final_data_object = eachDataObject
            else:
                final_data_object['run_metadata'] = pd.concat([final_data_object['run_metadata'], eachDataObject['run_metadata']]).reset_index(drop=True)
                if self.ignore_proteins is False:
                    final_data_object['protein_other_info'] = eachDataObject["protein_abundance"][["Accession"]]
                    final_data_object['protein_ID_Summary'] = pd.concat([final_data_object['protein_ID_Summary'], eachDataObject['protein_ID_Summary']]).reset_index(drop=True)
                if self.ignore_peptides is False:
                    final_data_object['peptide_other_info'] = eachDataObject["peptide_abundance"][["Annotated Sequence"]]
                    final_data_object['peptide_ID_Summary'] = pd.concat([final_data_object['peptide_ID_Summary'], eachDataObject['peptide_ID_Summary']]).reset_index(drop=True)
                duplicates_found = False
                
                #loop through to see if there are any duplicate files
                if self.ignore_proteins is False:
                    bad_cols = [eachCol for eachCol in final_data_object['protein_abundance'].loc[:, final_data_object['protein_abundance'].columns!='Accession'].columns if eachCol in eachDataObject["protein_abundance"]]
                    eachDataObject["protein_abundance"] =  eachDataObject["protein_abundance"].drop(bad_cols,axis=1)
                    
                    bad_cols = [eachCol for eachCol in final_data_object['protein_ID_matrix'].loc[:, final_data_object['protein_ID_matrix'].columns!='Accession'].columns if eachCol in eachDataObject["protein_ID_matrix"]]
                    eachDataObject["protein_ID_matrix"] =  eachDataObject["protein_ID_matrix"].drop(bad_cols,axis=1)
                    final_data_object['protein_abundance'] = pd.merge(final_data_object['protein_abundance'],eachDataObject['protein_abundance'],how="outer")
                    final_data_object['protein_ID_matrix'] = pd.merge(final_data_object['protein_ID_matrix'],eachDataObject['protein_ID_matrix'],how="outer")
                if self.ignore_peptides == False:
                    bad_cols = [eachCol for eachCol in final_data_object['peptide_abundance'].loc[:, final_data_object['peptide_abundance'].columns!='Annotated Sequence'].columns if eachCol in eachDataObject["peptide_abundance"]]
                    eachDataObject["peptide_abundance"] =  eachDataObject["peptide_ID_matrix"].drop(bad_cols,axis=1)
                    bad_cols = [eachCol for eachCol in final_data_object['peptide_ID_matrix'].loc[:, final_data_object['peptide_ID_matrix'].columns!='Annotated Sequence'].columns if eachCol in "peptide_ID_matrix"]
                    eachDataObject["peptide_ID_matrix"] =  eachDataObject["peptide_ID_matrix"].drop(bad_cols,axis=1)
                    final_data_object['peptide_abundance'] = pd.merge(final_data_object['peptide_abundance'],eachDataObject['peptide_abundance'],how="outer")
                    final_data_object['peptide_ID_matrix'] = pd.merge(final_data_object['peptide_ID_matrix'],eachDataObject['peptide_ID_matrix'],how="outer")
                    
        return final_data_object

    def sort_runs(self, data_object, settings_file):
        
        # Read the settings file
        settings_table = pd.read_table(settings_file, sep="\t")

        # Create a dictionary with the settings file {Condition: {column_head : value}}
        saved_settings = settings_table.set_index("Conditions").to_dict(orient="index")
        #any run with any of the filter_out items will not be used.

        # display(data_obj["protein_abundance"][pd.isna(data_obj["protein_abundance"])])
        # display(data_obj["protein_abundance"])
        # display(data_obj["protein_ID_Summary"])
        for eachGroup in saved_settings:
            i = 0

            # Create new column, "records"
            saved_settings[eachGroup]["records"] = []

            # Define what to filter out
            filterOutType = type(saved_settings[eachGroup]["filter_out"])
            if filterOutType == str or filterOutType == int or filterOutType == float and not pd.isna(saved_settings[eachGroup]["filter_out"]):
                filterOut = str.split(str(saved_settings[eachGroup]["filter_out"]),sep = ",")
            else:
                filterOut = ["M@di"]
            
            if len(str.split(str(saved_settings[eachGroup]["filter_in"]),sep = "@")) > 1: #multiple files, only some have the runs for this group
                user_list = {}
                #add all runs from all analyses to be probed
                for each_fileID in str.split(str(saved_settings[eachGroup]["filter_in"]),sep = "@")[1:]:
                    for eachIdentifier in data_object["run_metadata"]["Run Identifier"].drop_duplicates():    
                        currentRun = data_object["run_metadata"][data_object["run_metadata"]["Run Identifier"] == eachIdentifier]["Run Names"] 
                        if currentRun.size != 0:
                            if each_fileID == str.split(eachIdentifier,sep="-")[0] and list(currentRun)[0] not in user_list:
                                user_list[eachIdentifier] = list(currentRun)[0]
                        else:
                            print(currentRun)
                #filter for runs that relate to this gorup within those analyses
                # print(user_list)
                for run_id in user_list.keys():
                    run_name = user_list[run_id]
                    filter_ins = str.split(str.split(str(saved_settings[eachGroup]["filter_in"]),sep = "@")[0],",")
                    matches_all = True
                    for filter_in in filter_ins:
                        if bool(re.search(filter_in,run_name)) and (not any(item in run_name for item in filterOut)) and (run_name not in saved_settings[eachGroup]["records"]):
                            pass
                        else:
                            matches_all = False
                            break
                    if matches_all:
                        saved_settings[eachGroup]["records"].append(list(data_object["run_metadata"][data_object["run_metadata"]["Run Identifier"] == run_id]["Run Identifier"])[0]) 
                    else:
                        # print(filter_ins)
                        pass    
                    
            elif len(str.split(str(saved_settings[eachGroup]["filter_in"]),sep = "@")) == 1: #across all files or maybe there is only one
                    for run_name in data_object["run_metadata"]["Run Names"]:
                        run_id = data_object["run_metadata"].loc[data_object["run_metadata"]["Run Names"] ==  run_name,"Run Identifier"].tolist()[0]
                        filter_ins = str.split(str.split(str(saved_settings[eachGroup]["filter_in"]),sep = "@")[0],",")
                        matches_all = True
                        for filter_in in filter_ins:
                            if bool(re.search(filter_in,run_name)) and (not any(item in run_name for item in filterOut)) and (run_name not in saved_settings[eachGroup]["records"]):
                                pass
                            else:
                                matches_all = False
                                break
                        if matches_all:
                            saved_settings[eachGroup]["records"].append(list(data_object["run_metadata"][data_object["run_metadata"]["Run Identifier"] == run_id]["Run Identifier"])[0]) 
                        else:
                            # print(filter_ins)    
                            pass
                        
        #add the order of each column
        ignore_columns = ["filter_in","filter_out"]
        category_columns = [x for x in settings_table.columns.tolist() if x not in ignore_columns]

        for eachCol in category_columns:
            saved_settings["Order@"+eachCol] = settings_table[eachCol].drop_duplicates().tolist()
        
        return saved_settings

    def calculate_missing_values_MS2(self,data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Accession  3_TrypsinLysConly_3A4_channel2 3_TrypsinLysConly_3BC_channel1
    0     A0A096LP49                            0.00                                        10
    1     A0A0B4J2D5                        89850.26                                      3311
    2         A0AVT1                        83055.87                                    312312
        """
        if is_protein:
            name = "Accession"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        #initializes number of missing values to zero
        protein_columns = data_object[matrix_name].assign(missingValues=0)
        print(protein_columns)
        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name
        # print(data_object[matrix_name].columns)
        for each_column in data_object[matrix_name].loc[:, ~data_object[matrix_name].columns.str.contains(name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan,"NA":np.nan}) 
            protein_columns.loc[protein_columns[each_column] != "MS2", #ID/MBR are still missing values if you are only considering MS2
                                "missingValues"] += 1

            i += 1
        print(protein_columns)
        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        returnMatrix = pd.DataFrame({"Missing Values Rate": protein_columns["missingValuesRate"], name: protein_columns[name]})
        
        
        
        return returnMatrix

    def filter_by_missing_values(self,data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Accession  3_TrypsinLysConly_3A4_channel2 3_TrypsinLysConly_3BC_channel1
    0     A0A096LP49                            0.00                                        10
    1     A0A0B4J2D5                        89850.26                                      3311
    2         A0AVT1                        83055.87                                    312312
        """
        if is_protein:
            name = "Accession"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        #initializes a new column "missing values" with all rows equal to zero
        protein_columns = data_object[matrix_name].copy().assign(missingValues=0)

        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name

        # For each column associated with the group name (ex Helalib) except Accession/Annotated Sequence
        for each_column in data_object[matrix_name].loc[
                :, ~data_object[matrix_name].columns.str.contains(
                    name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan}) 

            #find missing values and increment those rows (a row is a protein/peptide) total number of missing values
            protein_columns.loc[(protein_columns[each_column] != "MS2")
                                &(protein_columns[each_column] != "MBR")
                                &(protein_columns[each_column] != "ID"), #this is more robust than using nan's in case something fails to convert
                                "missingValues"] += 1

            i += 1

        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        protein_columns = protein_columns.query(
            "missingValuesRate < @missing_value_thresh")
        
        data_object[abundance_name] = data_object[abundance_name].loc[data_object[abundance_name][name].isin (protein_columns[name].tolist())]
        data_object[matrix_name] = data_object[matrix_name].loc[data_object[matrix_name][name].isin (protein_columns[name].tolist())]
        data_object[other_info_name] = data_object[other_info_name].loc[data_object[other_info_name][name].isin( protein_columns[name].tolist())]

        # print(data_object[abundance_name])

        # protein_columns = protein_columns.loc[:,
        #                             protein_columns.columns.str.contains(name)]
        # protein_columns= protein_columns.dropna().drop_duplicates()
        # print(protein_columns)
        
        # # filter the data_object with the remaining proteins/peptides names
        # data_object[abundance_name] = protein_columns.merge(
        #     data_object[abundance_name],how='inner')
        # data_object[matrix_name] = protein_columns.merge(
        #     data_object[matrix_name], how = "inner")
        # data_object[other_info_name] = protein_columns.merge(
        #     data_object[other_info_name], how = "inner")
        # In case there is mismatch between ID table and abundance table,
        # mannually remove the row with all NaN values
        # keep rows in data_object[abundance_name] where at least two values are 
        # not NaN(do this to all rows except the first row), otherwise can't
        # calculate the stdev
        if ignore_nan:
            data_object[abundance_name] = data_object[abundance_name].dropna(
                thresh=2, subset=data_object[abundance_name].columns[1:])
            # This will cause the veen diagram to be different from R program
        
        return data_object

    def filter_by_missing_values_MS2(self,data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Accession  3_TrypsinLysConly_3A4_channel2
    0     A0A096LP49                            0.00
    1     A0A0B4J2D5                        89850.26
    2         A0AVT1                        83055.87
        """
        if is_protein:
            name = "Accession"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        protein_columns = data_object[matrix_name].copy()
        protein_columns["missingValues"] = 0
        # print(data_object[matrix_name].shape)
        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name
        for each_column in data_object[matrix_name].loc[:, ~data_object[matrix_name].columns.str.contains(name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan}) 
            protein_columns.loc[protein_columns[each_column] != "MS2", #ID/MBR are still missing values if you are only considering MS2
                                "missingValues"] += 1

            i += 1

        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        protein_columns = protein_columns.query(
            "missingValuesRate < @missing_value_thresh")
        
        data_object[abundance_name] = data_object[abundance_name].loc[data_object[abundance_name][name].isin (protein_columns[name].tolist())]
        data_object[matrix_name] = data_object[matrix_name].loc[data_object[matrix_name][name].isin (protein_columns[name].tolist())]
        data_object[other_info_name] = data_object[other_info_name].loc[data_object[other_info_name][name].isin( protein_columns[name].tolist())]

        # print(data_object[abundance_name])
        # protein_columns = protein_columns.loc[:,
        #                             protein_columns.columns.str.contains(name)]
        # protein_columns= protein_columns.dropna().drop_duplicates()
        # print(protein_columns)
        
        # # filter the data_object with the remaining proteins/peptides names
        # data_object[abundance_name] = protein_columns.merge(
        #     data_object[abundance_name], how = "inner")
        # data_object[matrix_name] = protein_columns.merge(
        #     data_object[matrix_name], how = "inner")
        # data_object[other_info_name] = protein_columns.merge(
        #     data_object[other_info_name], how = "inner")
        # # In case there is mismatch between ID table and abundance table,
        # # mannually remove the row with all NaN values
        # # keep rows in data_object[abundance_name] where at least two values are 
        # # not NaN(do this to all rows except the first row), otherwise can't
        # # calculate the stdev
        if ignore_nan:
            data_object[abundance_name] = data_object[abundance_name].dropna(
                thresh=2, subset=data_object[abundance_name].columns[1:])
            # This will cause the veen diagram to be different from R program
        # print(data_object[matrix_name])
        
        return data_object


    def NormalizeToMedian(self,abundance_data, apply_log2=False):
        """_Normalizes each column by multiplying each value in that column with
        the median of all values in abundances (all experiments) and then dividing
        by the median of that column (experiment)._
        we find applying log2 transform first gives more robust results for PCA etc.
        See https://pubs.acs.org/doi/10.1021/acsomega.0c02564
        Args:
            abundance_data (_pd_): _description_
            apply_log2 (_bool_,): _apply log2 to all result_.
        Returns:
            _type_: _description_
            format:
            'abundances':        Accession  3_TrypsinLysConly_3A4_channel2
            A0A096LP49                    0.000000e+00
        """
        # all the columns/sample list
        columns = [col for col in abundance_data.select_dtypes(include=[
                np.number])]
        data_matrix = abundance_data[columns].values
        # replace 0 with nan
        data_matrix[data_matrix == 0] = np.nan
        medianOfAll = np.nanmedian(data_matrix)
        
        #normalize all median, all median/current run all protein median
        # apply log2 to all the values if apply_log2 is True
        if apply_log2:    
            for each_column in columns:
                abundance_data[each_column] = (
                    np.log2(medianOfAll) * np.log2(abundance_data[each_column]) /
                    np.log2(np.nanmedian(abundance_data[
                        each_column].replace(0, np.nan))))
        else:
            for each_column in columns:
                abundance_data[each_column] = (
                    medianOfAll * abundance_data[each_column] /
                    np.nanmedian(abundance_data[
                        each_column].replace(0, np.nan)))
        #TODO divide by zero error encountered in log2, temporarily set to 0
        abundance_data = abundance_data.replace([np.inf, -np.inf], 0)

        return abundance_data

    def calculate_cvs(self,abundance_data):
        """_Calculate mean, stdev, cv for withn each protein/peptide abundance_

        Args:
            data_object (_type_): _full data frame_

        Returns:
            _type_: _df with Accession mean, stdev, cv for each protein/peptide_
        """
        if 'Accession' in abundance_data.columns:
            name = "Accession"
        if 'Annotated Sequence' in abundance_data.columns:
            name = "Annotated Sequence"
        abundance_data = abundance_data.assign(
            intensity=abundance_data.loc[:, ~abundance_data.columns.str.contains(
                name)].mean(axis=1, skipna=True),
            stdev=abundance_data.loc[:, ~abundance_data.columns.str.contains(
                name)].std(axis=1, skipna=True),
            CV=abs(abundance_data.loc[:, ~abundance_data.columns.str.contains(name)].std(
                axis=1, skipna=True) / abundance_data.loc[
                :, ~abundance_data.columns.str.contains(name)].mean(
                axis=1, skipna=True) * 100),
            count=abundance_data.loc[:, ~abundance_data.columns.str.contains(name)].count(axis=1))

        abundance_data = abundance_data.loc[:, [
                name, "intensity", "stdev", "CV","count"]]
        
        return abundance_data


    def t_test_from_summary_stats(self,m1, m2, n1, n2, s1, s2, equal_var=False):
        """_Calculate T-test from summary using ttest_ind_from_stats from
        scipy.stats package_

        Args:
            m1 (_type_): _mean list of sample 1_
            m2 (_type_): mean list of sample 2_
            n1 (_type_): sample size list of sample 1_
            n2 (_type_): sample size list of sample 2_
            s1 (_type_): standard deviation list of sample 1_
            s2 (_type_): standard deviation list of sample 2_
            equal_var (_type_, optional): False would perform Welch's
            t-test, while set it to True would perform Student's t-test. Defaults
            to False.

        Returns:
            _type_: _list of P values_
        """

        p_values = []
        for i in range(len(m1)):
            _, benjamini = ttest_ind_from_stats(
                m1[i], s1[i], n1[i], m2[i], s2[i], n2[i], equal_var=equal_var)
            p_values.append(benjamini)

        return p_values
    
    def percentile(SELF,percentile):
        def percentile_(x):
            return x.quantile(percentile)
        percentile_.__name__ = 'percentile_{:02.0f}'.format(percentile*100)
        return percentile_

    def impute_knn(self,abundance_data, k=5):
        """_inpute missing value from neighbor values_

        Args:
            abundance_data (_type_): _description_
            k (int, optional): _number of neighbors used_. Defaults to 5.
        Returns:
            _type_: _description_
            TODO: this knn imputer produces slightly different results (about 4%)
            from the one in R. Need to figure out why
        """
        name = abundance_data.columns[0]

        names = abundance_data[name]
        # x = abundance_data.select_dtypes(include=['float64', 'int64'])
        # imputer = KNNImputer(n_neighbors=k)
        # x_imputed = pd.DataFrame(imputer.fit_transform(x), columns=x.columns)


        x = abundance_data.select_dtypes(include=['float', 'int'])
        # print(abundance_data)

        imputer = KNNImputer(n_neighbors=k)
        x_imputed = imputer.fit_transform(x)
        x_imputed = pd.DataFrame(x_imputed, columns=x.columns)
        # print(x_imputed.columns)

        # Replace the original values in abundance_data with imputed values
        # x_imputed.insert(loc=0,column=name,value=names)

        # Check if the column already exists
        if name in x_imputed.columns:
            # Update the existing column
            x_imputed[name] = names
        else:
            # Insert the new column
            x_imputed.insert(loc=0, column=name, value=names)
        
        return x_imputed

    def CalculatePCA(self,abundance_object, infotib,log2T = False):
        """_inpute PCA transformed and variance explained by each principal
        component_
        """
        name = abundance_object.columns[0]
        x = abundance_object
        
        sampleNames = x.columns[~x.columns.str.contains(
            name)].to_frame(index=False)

        if log2T: #apply log2 transformation
            x = np.log2(x.loc[:, ~x.columns.str.contains(name)].T.values)
        else:
            x = x.loc[:, ~x.columns.str.contains(name)].T.values
        # filter out columns with all zeros
        is_finite_col = np.isfinite(np.sum(x, axis=0))
        x_filtered = x[:, is_finite_col]

        
        # Instantiate PCA   
        pca = PCA()
        #
        # Determine transformed features
        #
        x_pca = pca.fit_transform(x_filtered)
        #
        # Determine explained variance using explained_variance_ration_ attribute
        #
        exp_var_pca = pca.explained_variance_ratio_
        #
        # Cumulative sum of eigenvalues; This will be used to create step plot
        # for visualizing the variance explained by each principal component.
        #
        cum_sum_eigenvalues = np.cumsum(exp_var_pca)
        #
        # convert numpy array to pandas dataframe for plotting
        
        pca_panda = pd.DataFrame(x_pca, columns=[
            'PC' + str(i+1) for i in range(x_pca.shape[1])])
        # add sample names to the dataframe
        pca_panda = pd.concat(
            [infotib, pca_panda], axis=1, join='inner')
        
        return pca_panda, exp_var_pca


    
# def filter_by_name(self,data_dict, runname_list):
#         """_Filter the data_dict based on runname_list, only keep the columns
#         of the data_dict that are in the runname_list_
#         Args:

#         Returns:
#             _type_: _description_
#         """

#         # make dict for each runname, no Accession/sequence
#         nameDict_channels = dict(zip(data_dict["run_metadata"]["Run Names"],data_dict["run_metadata"]["Channel Identifier"]))

#         nameDict_runs = dict(zip(data_dict["run_metadata"]["Run Names"],data_dict["run_metadata"]["Run Identifier"]))
        
#         identifier_list = []
        
#         identifier_list_plus = []

#         run_id_list = []

#         if "Annotated Sequence" in runname_list:
#             runname_list.remove("Annotated Sequence")
#         if "Accession" in runname_list:
#             runname_list.remove("Accession")
#         for eachName in runname_list:
#             run_id_list.append(nameDict_runs[eachName])
#         for eachName in runname_list:
#             identifier_list.append(nameDict_channels[eachName])
#         for eachName in runname_list:
#             identifier_list_plus.append(nameDict_channels[eachName])


#         filtered_data = {}
#     # filtered_data["meta"] = data_dict["meta"]
#         runname_list.extend(["Annotated Sequence","Accession"])
#         identifier_list_plus.extend(["Annotated Sequence","Accession"])

#         #filtered_data["run_metadata"] = [item for item in data_dict[
#         #   "run_metadata"] if item in runname_list]
        
#         filtered_data["run_metadata"] = data_dict["run_metadata"][
#             data_dict["run_metadata"]["Run Names"].isin(
#                 runname_list)]  
#         filtered_data["protein_abundance"] = data_dict["protein_abundance"][[
#             col for col in data_dict["protein_abundance"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["peptide_abundance"] = data_dict["peptide_abundance"][[
#             col for col in data_dict["peptide_abundance"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["protein_other_info"] = data_dict["protein_other_info"][[
#             col for col in data_dict["protein_other_info"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["peptide_other_info"] = data_dict["peptide_other_info"][[
#             col for col in data_dict["peptide_other_info"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["protein_ID_matrix"] = data_dict["protein_ID_matrix"][[
#             col for col in data_dict["protein_ID_matrix"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["peptide_ID_matrix"] = data_dict["peptide_ID_matrix"][[
#             col for col in data_dict["peptide_ID_matrix"].columns if any(
#                 word == col for word in identifier_list_plus)]]
#         filtered_data["protein_ID_Summary"] = data_dict["protein_ID_Summary"][
#             data_dict["protein_ID_Summary"]["names"].isin(
#                 run_id_list)]
#         filtered_data["peptide_ID_Summary"] = data_dict["peptide_ID_Summary"][
#             data_dict["peptide_ID_Summary"]["names"].isin(
#                 run_id_list)]
#         return filtered_data

    def filter_by_id(self,data_dict, run_id_list):
        """_Filter the data_dict based on runname_list, only keep the columns
        of the data_dict that are in the runname_list_
        Args:

        Returns:
            _type_: _description_
        """

        # make dict for each runname, no Accession/sequence
        
        nameDict_channels = dict(zip(data_dict["run_metadata"]["Run Identifier"],data_dict["run_metadata"]["Run Identifier"]))

        identifier_list = []
        
        identifier_list_plus = []


        if "Annotated Sequence" in run_id_list:
            run_id_list.remove("Annotated Sequence")
        if "Accession" in run_id_list:
            run_id_list.remove("Accession")
        for eachName in run_id_list:
            identifier_list.append(nameDict_channels[eachName])
        for eachName in run_id_list:
            identifier_list_plus.append(nameDict_channels[eachName])


        filtered_data = {}
        # filtered_data["meta"] = data_dict["meta"]
        run_id_list.extend(["Annotated Sequence","Accession"])
        identifier_list_plus.extend(["Annotated Sequence","Accession"])

        #filtered_data["run_metadata"] = [item for item in data_dict[
        #   "run_metadata"] if item in runname_list]
        
        filtered_data["run_metadata"] = data_dict["run_metadata"][
            data_dict["run_metadata"]["Run Identifier"].isin(
                run_id_list)]  
        if self.ignore_peptides == False:
            filtered_data["peptide_abundance"] = data_dict["peptide_abundance"][[
                col for col in data_dict["peptide_abundance"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_other_info"] = data_dict["peptide_other_info"][[
                col for col in data_dict["peptide_other_info"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_ID_matrix"] = data_dict["peptide_ID_matrix"][[
                col for col in data_dict["peptide_ID_matrix"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_ID_Summary"] = data_dict["peptide_ID_Summary"][
                data_dict["peptide_ID_Summary"]["names"].isin(
                    run_id_list)]
        if self.ignore_proteins == False:
            filtered_data["protein_abundance"] = data_dict["protein_abundance"][[
                col for col in data_dict["protein_abundance"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_other_info"] = data_dict["protein_other_info"][[
                col for col in data_dict["protein_other_info"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_ID_matrix"] = data_dict["protein_ID_matrix"][[
                col for col in data_dict["protein_ID_matrix"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_ID_Summary"] = data_dict["protein_ID_Summary"][
                data_dict["protein_ID_Summary"]["names"].isin(
                    run_id_list)]
      
        return filtered_data
        
    def filter_by_channel_id(self,data_dict, run_id_list):
        """_Filter the data_dict based on runname_list, only keep the columns
        of the data_dict that are in the runname_list_
        Args:

        Returns:
            _type_: _description_
        """

        # make dict for each runname, no Accession/sequence
        if self.data_type == "TMT":
            nameDict_channels = dict(zip(data_dict["run_metadata"]["Channel Identifier"],data_dict["run_metadata"]["Channel Identifier"]))
        elif self.data_type == "LF":
            nameDict_channels = dict(zip(data_dict["run_metadata"]["Run Identifier"],data_dict["run_metadata"]["Run Identifier"]))

        identifier_list = []
        
        identifier_list_plus = []


        if "Annotated Sequence" in run_id_list:
            run_id_list.remove("Annotated Sequence")
        if "Accession" in run_id_list:
            run_id_list.remove("Accession")
        for eachName in run_id_list:
            identifier_list.append(nameDict_channels[eachName])
        for eachName in run_id_list:
            identifier_list_plus.append(nameDict_channels[eachName])


        filtered_data = {}
    # filtered_data["meta"] = data_dict["meta"]
        run_id_list.extend(["Annotated Sequence","Accession"])
        identifier_list_plus.extend(["Annotated Sequence","Accession"])

        #filtered_data["run_metadata"] = [item for item in data_dict[
        #   "run_metadata"] if item in runname_list]
        
        filtered_data["run_metadata"] = data_dict["run_metadata"][
            data_dict["run_metadata"]["Run Identifier"].isin(
                run_id_list)]  
        if self.ignore_peptides == False:
            filtered_data["peptide_abundance"] = data_dict["peptide_abundance"][[
                col for col in data_dict["peptide_abundance"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_other_info"] = data_dict["peptide_other_info"][[
                col for col in data_dict["peptide_other_info"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_ID_matrix"] = data_dict["peptide_ID_matrix"][[
                col for col in data_dict["peptide_ID_matrix"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["peptide_ID_Summary"] = data_dict["peptide_ID_Summary"][
                data_dict["peptide_ID_Summary"]["names"].isin(
                    run_id_list)]
        if self.ignore_proteins == False:
            filtered_data["protein_abundance"] = data_dict["protein_abundance"][[
                col for col in data_dict["protein_abundance"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_other_info"] = data_dict["protein_other_info"][[
                col for col in data_dict["protein_other_info"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_ID_matrix"] = data_dict["protein_ID_matrix"][[
                col for col in data_dict["protein_ID_matrix"].columns if any(
                    word == col for word in identifier_list_plus)]]
            filtered_data["protein_ID_Summary"] = data_dict["protein_ID_Summary"][
                data_dict["protein_ID_Summary"]["names"].isin(
                    run_id_list)]
        return filtered_data
        
