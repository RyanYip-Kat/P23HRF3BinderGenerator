import json
import os
import pandas as pd
import shutil
from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
from tqdm import tqdm

import warnings
# 全局忽略所有警告
warnings.filterwarnings('ignore')

def cif_to_pdb_biopython(cif_file, pdb_file):
    # 使用 MMCIFParser 读取 CIF
    parser = PDB.MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure('struct', cif_file)

        # 使用 PDBIO 保存为 PDB
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        print(f"Biopython 转换成功: {pdb_file}")
    except Exception as e:
        print(f"Biopython 转换失败: {e}")

def get_sequence_from_pdb_chain(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    model = structure[0]
    chain_names =[x.id for x in model.get_chains()]
    ppb = PDB.PPBuilder()
    result = {chain:"" for chain in chain_names}
    for chain_name in chain_names:
        pp =ppb.build_peptides(model[chain_name])
        for pp in pp:
            result[chain_name]+=str(pp.get_sequence())
    return result

def get_sequence_from_cif_chain(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure('CIF_structure', cif_file)
    model = structure[0]
    chain_names =[x.id for x in model.get_chains()]
    ppb = PDB.PPBuilder()
    result = {chain:"" for chain in chain_names}
    for chain_name in chain_names:
        pp =ppb.build_peptides(model[chain_name])
        for pp in pp:
            result[chain_name]+=str(pp.get_sequence())
    return result


batch_id = "batch1"
original_pdb_file = "/home/data1/ryanyip/project/RFdiffusion3/liuzhong/P23H_opsin_chainA.pdb"
result_dir = f"liuzhong/protein_binder/P23H_V2/{batch_id}"
result_csv_file = f'{result_dir}/results.csv'
seq_ranges = [70,86]
prefix=f"{batch_id}_p23h"
results =[]

diff_generate_refold = 0

target_seq = list(get_sequence_from_pdb_chain(original_pdb_file).values())[0]
for seq_len  in range(*seq_ranges):
    out_json_file = f'{result_dir}/{seq_len}/pipeline_results.json'
    if not os.path.exists(out_json_file):
        continue
    out_dir = os.path.dirname(out_json_file)
    dataset = json.load(open(out_json_file,encoding="utf-8"))
    print(f'{len(dataset)} results')

    outs =[]
    for data in tqdm(dataset, desc='Processing results',total=len(dataset)):
        struct_idx = data['structure_id']
        mpnn_section_list = data['section2']
        mpnn_output_list=[x['output_dict'] for x in mpnn_section_list['designed_sequences']]
        mpnn_sequence = [x['sequence'] for x in mpnn_section_list['designed_sequences']]
        rf3_section_list = data['section3']['rf3_predictions']
        
        final_section_list = data['section4']
        validation_results = final_section_list['validation_results']

        # print(f'ProteinMPNN 优化了 {len(mpnn_sequence)} 个序列')
        assert len(mpnn_output_list) == len(rf3_section_list) == len(validation_results)

        cif_dir = os.path.join(out_dir, struct_idx)
        for idx in range(len(mpnn_output_list)):
            rfd3_cif_file = os.path.join(cif_dir, f'sequence_{idx}',"generated.cif")
            rf3_cif_file = os.path.join(cif_dir, f'sequence_{idx}',"refolded.cif")
            rf3_metric = rf3_section_list[idx]['best_model_summary']
            section4_metric = validation_results[idx]

            mpnn_sequence_recovery = mpnn_output_list[idx]['sequence_recovery']
                
            rfd3_sequence_dict = get_sequence_from_cif_chain(rfd3_cif_file)
            rf3_sequence_dict = get_sequence_from_cif_chain(rf3_cif_file)

            # assert rfd3_sequence_dict['B'] == rf3_sequence_dict['B']
            if rfd3_sequence_dict['B'] != rf3_sequence_dict['B']:
                diff_generate_refold +=1
                continue
            # target_seq = rfd3_sequence_dict['B']
            out = {
                'peptide_len': seq_len,
                'structure_id': struct_idx,
                'mpnn_id': idx,
                'mpnn_sequence_recovery': mpnn_sequence_recovery,
                'final_sequence': mpnn_sequence[idx][:-len(target_seq)],
                'rfd3_cif_file': rfd3_cif_file,
                'rfd3_binder_seq': rfd3_sequence_dict['A'],
                'rf3_cif_file': rf3_cif_file,
                'rf3_binder_seq': rf3_sequence_dict['A'],
                }
            if rf3_metric is not None:
                for k,v in rf3_metric.items():
                    out[k] = v
            if section4_metric is not None:
                for k in ['backbone_rmsd', 'interpretation', 'comparison_atoms']:
                    out[k] = section4_metric[k]
            outs.append(out)
    results.extend(outs)

df = pd.DataFrame(results)
# df['target_seq'] = target_seq
df.to_csv(result_csv_file,index=False)
print(df.head())
print(f'Notice: {len(df)} results saved to {result_csv_file}')
print(f'diff_generate_refold: {diff_generate_refold}')

#  copy file
final_csv_file = f'{result_dir}/results_final.csv'
output_dir = f"{result_dir}/processed_results"
rfd3_output_dir = f"{output_dir}/rfd3"
rf3_output_dir = f"{output_dir}/rf3"
os.makedirs(rf3_output_dir, exist_ok=True)
os.makedirs(rfd3_output_dir, exist_ok=True)

dataset=pd.read_csv(result_csv_file)
outs =[]
for index,row in tqdm(dataset.iterrows(),total=len(dataset),desc="Processing"):
    rfd3_cif_file = row['rfd3_cif_file'] # generated.cif
    rf3_cif_file = row['rf3_cif_file'] # refolded.cif
    peptide_len = int(row['peptide_len'])

    mpnn_id = str(row['mpnn_id'])
    structure_id = str(row['structure_id'])

    # rfd3_dst_file = f"{rfd3_output_dir}/{prefix}_{structure_id}_{mpnn_id}_{peptide_len}_generated.cif"
    # rf3_dst_file = f"{rf3_output_dir}/{prefix}_{structure_id}_{mpnn_id}_{peptide_len}_refolded.cif"
    rfd3_dst_file = f"{rfd3_output_dir}/{prefix}_{structure_id}_{mpnn_id}_{peptide_len}_generated.pdb"
    rf3_dst_file = f"{rf3_output_dir}/{prefix}_{structure_id}_{mpnn_id}_{peptide_len}_refolded.pdb"
    
    row['rfd3_file'] = rfd3_dst_file
    row['rf3_file'] = rf3_dst_file
    outs.append(row.to_dict())

    # shutil.copyfile(rfd3_cif_file, rfd3_dst_file)
    # shutil.copyfile(rf3_cif_file, rf3_dst_file)
    cif_to_pdb_biopython(rfd3_cif_file, rfd3_dst_file)
    cif_to_pdb_biopython(rf3_cif_file, rf3_dst_file)

df = pd.DataFrame(outs)
df.to_csv(final_csv_file,index=False)
print("Done")
print(df.head())
