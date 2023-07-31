# These are the functions used to clean the Uniprot database

# Fetch all uniprot data
import requests
import pandas as pd
import numpy as np 

class UniprotDataFetcher:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id 
    
    def get_sequence_length(self):
        url = f"https://www.uniprot.org/uniprot/{self.uniprot_id}.fasta"
        response = requests.get(url)
        if response.ok:
            sequence = response.text.split('\n', 1)[1].replace('\n', '')
            return len(sequence)
        else:
            print("Failed to retrieve sequence.")
            return None
    
    def get_data(self):
        # url for specific unirpot ID data
        url = f'https://rest.uniprot.org/uniprot/{self.uniprot_id}'
        # Fetch json response
        data = requests.get(url).json()
        # only get uniprot data
        try:
            d = data["uniProtKBCrossReferences"]
        except:
            raise ValueError("Could not get Uniprot data...")
        return d 
    # Given all the uniprot data of a specifc uniprot ID 
    # return information concerning the PDB files associated with the ID
    def get_pdb_info(self):
        d = self.get_data()
        pdb_info = [] 
        for item in d:
            if item["database"] == "PDB":
                pdb_info.append(item) 
        return pdb_info
    
    def extract_chain_identifiers(self, chains):
        chain_identifiers = [] 
        # split the input string by comma 
        pairs = chains.split(",")
        
        for pair in pairs:
            # remove leading/trailing spaces
            pair = pair.strip() 
            # Extract the cahin identifer
            chain_id = pair[0] 
            # append id to chain identifiers list
            chain_identifiers.append(chain_id)
            
        return chain_identifiers
            
    def remove_rows_with_nmr(self, df):
        df_filtered = df.drop(df[df["Method"] == "NMR"].index)
        # If NMR structures are the only structures for the specific Uniprot ID then keep original df
        if len(df_filtered) == 0:
            return df
        return df_filtered
    
    def create_empty_df(self):
        data = {
            'PDB ID': [self.uniprot_id],
            'Resolution': [np.nan],
            'Locations' : locations,
            'Chain': [np.nan],
            'PDB_Sequence_Length': [np.nan],
            #'UniProt_Sequence_Length': [np.nan],
            'Method': [np.nan]
        }
        df = pd.DataFrame(data)
        return df
    
    def get_chain_only(self, pdb):
        structure = self.get_pdb_info()
        chain = ""
        for item in structure:
            if item['database'] =='PDB':
                if item['id'] == pdb:
                    properties = item['properties']
                    for prop in properties:
                        if prop['key'] == "Chains":
                            chain = prop["value"]
        print(f"Chain:{chain}")
        return chain
            
    def convert_structure_to_dataframe(self):
        structure = self.get_pdb_info()
        # Edge Case: No PDB files for Uniprot ID
        if len(structure) == 0:
            empty_df = self.create_empty_df()
            return empty_df
            
        # Initialize empty lists to store the extracted values
        pdb_ids = []
        resolutions = []
        chains = []
        sequence_lengths = []
        methods = []
        locations = []

        for item in structure:
            if item['database'] == 'PDB':
                pdb_id = item['id']
                # Edge Case - if there are no pdb files for the specific Uniprot ID 
                properties = item['properties']
                for prop in properties:
                    if prop['key'] == 'Resolution':
                        # Extract only the number from the resolution
                        resolution = prop["value"].split()[0]
                        if resolution != '-':
                            resolution = float(resolution)
                    elif prop['key'] == 'Chains':
                        # set coverage to 0
                        cov = 0
                        chain = prop['value']
                        locations.append(chain)
                        # indicating multiple different domains with different length corespondong to molecule
                        
                        if ',' in chain:
                            # coverage value
                            values = [item.split('=')[1] for item in chain.split(',')]
                            for val in values: 
                                start, end = map(int, val.split('-'))
                                cov += end - start + 1
                        else: 
                            try:
                                start, end = map(int, chain.split("=")[1].split('-'))
                                cov = end - start + 1
                            except:
                                cov = 0
                            
                    else:
                        # get method of extraction
                        method = prop["value"]

                # Append the extracted values to their respective lists
                pdb_ids.append(pdb_id)
                resolutions.append(resolution)
                chains.append(','.join(self.extract_chain_identifiers(chain)))
                sequence_lengths.append(cov)
                methods.append(method)
                #uniprot_sequence_length = [self.get_sequence_length() for _ in range(len(pdb_ids))]
                

        # Create a DataFrame from the extracted data
        data = {
            'PDB ID': pdb_ids,
            'Resolution': resolutions,
            'Locations' : locations,
            'Chain': chains,
            'PDB_Sequence_Length': sequence_lengths,
            #'UniProt_Sequence_Length' : uniprot_sequence_length,
            "Method": methods
        }

        df = pd.DataFrame(data)
        df = self.remove_rows_with_nmr(df)
        
        # Remove any dataframes where Chain is a numberic value (indicating AI generated structures??) 
        df_c = df.copy() 
        df_c["Chain"] = pd.to_numeric(df_c["Chain"], errors ="coerce")
        df_filtered = df[df_c["Chain"].isna()] 
        df = df_filtered
        return df
    
                
class PDB_analyzer:
    def __init__(self, pdb_id, chains):
        self.pdb_id = pdb_id
        self.chains = chains
    
    def fetch_pdb_content(self):
        url = f'https://files.rcsb.org/download/{self.pdb_id}.pdb'
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            return np.nan
    
    def does_PDB_FILE_Exist(self):
        if isinstance(self.fetch_pdb_content(),str):
            return True
        else:
            return False
    
    def adp_file_path(self):
        if self.fetch_pdb_content() == np.nan:
            return np.nan
        atom_file_path = f"pdb_files/{self.pdb_id}.apdb"
        return atom_file_path 
    
    def write_filter_alpha_carbons(self):
        if self.adp_file_path() == np.nan:
            return np.nan
        # initialize both atom count and alpha carbon count to compare later
        atom_count = 0
        ca_count = 0 
        # Initialize where we will store atom files
        apdb_file_path = self.adp_file_path()
        pdb_content = self.fetch_pdb_content()
        #print(f"FILE: {apdb_file_path}, Length of pdb content {len(pdb_content)}")
        if pdb_content == np.nan: 
            return np.nan
        with open(apdb_file_path, 'w') as apdb_file:
           # print(pdb_content)
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM'):
                    atom_name = line[12:16].strip()
                    apdb_file.write(line + '\n')
                    atom_count += 1
                    if atom_name == 'CA':
                        ca_count += 1
        
        if ca_count == atom_count:
            return True
        else:
            return False
        
    # Given the chain ID, the PDB file: get the coverage of the pdb file
    def get_chain_coverage_depth(self, chain):
        atom_file_path = self.adp_file_path()
        if atom_file_path == np.nan:
            return np.nan, np.nan
        chain_id = chain
        first_residue_number = float('inf') 
        last_residue_number = float('-inf')
        count = 0
        curr = 0
        try: 
            with open(atom_file_path, 'r') as atom_file:
                
                for line in atom_file:
                    if line.startswith('ATOM'):
                        current_chain_id = line[21]
                        if current_chain_id == chain_id:
                            residue_number = int(line[22:26].strip())
                            if curr != residue_number:
                                curr = residue_number
                                count += 1 
                                # Debugging
                                #print(f"Current: {curr}")
                                #print(f"Residue Numbe: {residue_number}")
                                #print(f"Count: {count}")
                                #print("------------------------------------------------------------")
                            first_residue_number = min(first_residue_number, residue_number)
                            last_residue_number = max(last_residue_number, residue_number)
                        #print(first_residue_number)
                        #print(last_residue_number)
        except: 
            return np.nan, np.nan
                        

        if first_residue_number == float('inf') or last_residue_number == float('-inf'):
            #print(f"Count = {count} ")
            #print(f"Atom File Path = {atom_file_path}")
            #print(f"Chain ID: {self.chains} ") 
            return np.nan, np.nan
            #raise Exception(f"No ATOM lines found for chain {atom_file_path} in the file.")

        coverage = count - 1
        depth = last_residue_number - first_residue_number

        return coverage, depth
    
    def get_coverage_and_CA_info(self):
        #print(f"{self.pdb_id} : does it exist : {self.does_PDB_FILE_Exist()}")
        x = self.does_PDB_FILE_Exist()
        if x == False:
            return np.nan, np.nan, np.nan
        else:
            # get whether the file contains all alpha carbons or not 
            check = self.write_filter_alpha_carbons()
            # split the chain based on comma 
            chains = ','.join(self.chains)
            chains = chains.split(',')
            chains = list(filter(lambda value: value != '', chains))
            #print(f"Final function chains: {chains}")
            coverage = 0
            depth = 0 
            for chain in chains: 
                cov, dep = self.get_chain_coverage_depth(chain)
                coverage += cov 
                depth += dep
            return coverage, depth, check
    
         
        

    
            
               
        
    
    
    

    
            
        