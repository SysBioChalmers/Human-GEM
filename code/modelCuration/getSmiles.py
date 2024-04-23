import pandas as pd
df = pd.read_excel('./Human-GEM.xlsx',sheet_name='METS')
df1 = pd.read_excel('./metabolites_smiles_.xlsx')
df2 = pd.merge(df,df1,on='REPLACEMENT ID',how='left')
df2.to_excel('./metabolites_smiles.xlsx',index=False)


import pandas as pd
df = pd.read_excel('./metabolites_smiles.xlsx')
df1 = pd.read_csv('./database/mnx_chem_depr.tsv', sep='\t')
id = df1['deprecated_ID'].tolist() 
df['metMetaNetXID_new'] = ''
for i in range(len(df)):
    if df['metMetaNetXID'][i] in id:
        index = id.index(df['metMetaNetXID'][i])
        df['metMetaNetXID_new'][i] = df1['ID'][index]
df.to_excel('./metabolites_smiles.xlsx', index=False)


#get SMILES from database and model
import pandas as pd
from tqdm import tqdm
#df = pd.read_csv('./database/chebi_id_smiles.csv')
df = pd.read_csv('./database/kegg_compound.txt', sep='\t')
#df = pd.read_csv('./database/chebi_second_id_smiles.csv')
#df = pd.read_csv('./database/recon3d_smiles.csv')

#id = df['ChEBI ID'].to_list()
id = df['KEGG'].to_list()
#id = df['Secondary ChEBI ID'].to_list()
#id = df['metRecon3DID'].to_list()
df1 = pd.read_excel('./metabolites_smiles.xlsx')
#s = 'metChEBIID'
s = 'metKEGGID'
#s = 'metRecon3DID'
#df1['SMILES'] = None
for i in tqdm(range(len(df1)),total=len(df1)):
    if df1['SMILES'].isna()[i] == False:
        continue
    try:
        if i == 0:
            if df1[s].isna()[i] == True:
                continue
            else:
                if df1[s][i] in id:
                    index = id.index(df1[s][i])
                    df1['SMILES'][i] = df['SMILES'][index]
                
        else:
            if df1[s][i] == df1[s][i-1]:
                continue
            else:
                if df1[s].isna()[i] == True:
                    continue
                else:
                    if df1[s][i] in id:
                        index = id.index(df1[s][i])
                        df1['SMILES'][i] = df['SMILES'][index]
    except:
        continue
      
for i in tqdm(range(len(df1)),total=len(df1)):
    if i != 0:
        if df1[s][i] == df1[s][i-1]:
            df1['SMILES'][i] = df1['SMILES'][i-1]

df1.to_excel('metabolites_smiles.xlsx', index=False)


#get SMILES from database and model
import pandas as pd
from tqdm import tqdm

df = pd.read_csv('./database/mnx_chem_prop.tsv', sep='\t')
id = df['id'].to_list()

df1 = pd.read_excel('./metabolites_smiles.xlsx')
s = 'metMetaNetXID'
for i in tqdm(range(len(df1)),total=len(df1)):
    if df1['SMILES'].isna()[i] == False:
        continue
    try:
        if i == 0:
            if df1[s].isna()[i] == True:
                continue
            else:
                for j in df1[s][i].split(';'):
                    if j in id:
                        index = id.index(j)
                        df1['SMILES'][i] = df['SMILES'][index]
                        break
                
        else:
            if df1[s][i] == df1[s][i-1]:
                continue
            else:
                if df1[s].isna()[i] == True:
                    continue
                else:
                
                    for j in df1[s][i].split(';'):
                        if j in id:
                            index = id.index(j)
                            df1['SMILES'][i] = df['SMILES'][index]
                            break
    except:
        continue
      
for i in tqdm(range(len(df1)),total=len(df1)):
    if i != 0:
        if df1[s][i] == df1[s][i-1]:
            df1['SMILES'][i] = df1['SMILES'][i-1]

df1.to_excel('metabolites_smiles.xlsx', index=False)


import pubchempy as pcp

import pandas as pd
from tqdm import tqdm
failed = []
df1 = pd.read_excel('./metabolites_smiles.xlsx')
#df1['SMILES'] = ''
for i in tqdm(range(len(df1)),total=len(df1)):
    if df1['SMILES'].isna()[i] == False:
        continue
    try:
        #df1['metChEBIID'][i] = df1['metChEBIID'][i].replace('CHEBI:','')
        if i == 0:
            if df1['metPubChemID'].isna()[i] == True:
                continue
            else:
                cid = int(df1['metPubChemID'][i])
                df1['SMILES'][i] = pcp.Compound.from_cid(cid).isomeric_smiles
        else:
            if df1['metsNoComp'][i] == df1['metsNoComp'][i-1]:
                continue
            else:
                if df1['metPubChemID'].isna()[i] == True:
                    continue
                else:
                    cid = int(df1['metPubChemID'][i])
                    df1['SMILES'][i] = pcp.Compound.from_cid(cid).isomeric_smiles
                    #print(df1['SMILES'][i])
    except:
        failed.append(df1['metPubChemID'][i])
        continue
for i in tqdm(range(len(df1)),total=len(df1)):
    if i != 0:
        if df1[s][i] == df1[s][i-1]:
            df1['SMILES'][i] = df1['SMILES'][i-1]
df1.to_excel('./metabolites_smiles.xlsx', index=False)
# Get the SMILES string of the compound


from rdkit import Chem
def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)
df = pd.read_excel('metabolites_smiles.xlsx')
df['standard_smiles'] = df['SMILES'].apply(standardize_smiles)
df.to_excel('metabolites_smiles.xlsx', index=False)



filepath = 'D:\\All_Human_GTEx\\'
df_model_mets = pd.read_excel(filepath+'metabolites.xlsx')
# df_HMDB_mets = pd.read_excel(filepath+'副本hmdb_metabolites_217920_orign.xlsx')

mets_id1 = df_model_mets['mets'].values.tolist()
mets_id2 = df_model_mets['metsNoComp'].values.tolist()
mets_SMILE = df_model_mets['SMILES'].values.tolist()

mets_inchikey = []
mets_inchi = []
for i in mets_SMILE:
    i = str(i)
    if i == 'nan':
        mets_inchikey.append('')
        mets_inchi.append('')

    elif i.startswith('*'):
        mets_inchikey.append('')
        mets_inchi.append('')
        
    else:
        # print(i)
        mol = Chem.MolFromSmiles(i)
        mets_inchikey.append(Chem.MolToInchiKey(mol))
        mets_inchi.append(Chem.MolToInchi(mol))

filepath1 = 'D:\\All_Human_GTEx\\'
output = {'mets':mets_id1, 'metsNoComp':mets_id2, 'SMILES': mets_SMILE, 'inchikey':mets_inchikey, 'inchi':mets_inchi}
output_tsv = pd.DataFrame(output)
output_tsv.to_csv(filepath1+'metabolites_SMILES_Inchi.tsv', sep="\t", index=False)
