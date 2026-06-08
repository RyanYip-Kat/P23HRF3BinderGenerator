import warnings
import os
import time
import numpy as np
import json
from datetime import datetime
warnings.filterwarnings('ignore', module='atomworks')
warnings.filterwarnings('ignore')

from atomworks.io.utils.io_utils import to_cif_file

from lightning.fabric import seed_everything
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from biotite.structure import get_residue_starts
from biotite.sequence import ProteinSequence


from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput

from biotite.structure import rmsd, superimpose
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES

from tqdm import tqdm


# =========================================================
# 设定推理参数
seed = 3872541 # Set seed for reproducibility
length_min = 70 # Set minimum design length
length_max = 86 # Set maximum design length
diffusion_batch_size=10 # Generate "diffusion_batch_size" structures per batch
rfd3_n_batches = 5 # Generate "rfd3_n_batches" batches to run RFD3
mpnn_batch_size=10 # Generate "mpnn_batch_size" structures for mpnn

# 总生成的结构数量 ： diffusion_batch_size * rfd3_n_batches * mpnn_batch_size * （length_max - length_min）
total_structures = diffusion_batch_size * rfd3_n_batches * mpnn_batch_size * (length_max - length_min)
print(f'------ 总共要生成{total_structures}个结构 ------')

start = time.time()
# =========================================================

input_pdb = "/home/data1/ryanyip/project/RFdiffusion3/liuzhong/P23H_opsin_chainA.pdb"
global_prefix = 'p23h'
out_path="liuzhong/protein_binder/P23H_V2/batch1" 
os.makedirs(out_path,exist_ok=True)

# =========================================================
rfd3_ckpt_dir = "checkpoint/rfd3_latest.ckpt"
mpnn_ckpt_dir = "checkpoint/ligandmpnn_v_32_010_25.pt"
rf3_ckpt_dir = "checkpoint/rf3_foundry_01_24_latest_remapped.ckpt"


# See mpnn.utils.inference.MPNN_GLOBAL_INFERENCE_DEFAULTS for all options
engine_config = {
    "model_type": "ligand_mpnn",  # or "protein_mpnn" for vanilla ProteinMPNN
    "is_legacy_weights": True,    # Required for now for ligand_mpnn and protein_mpnn
    "out_directory": None,        # Return results in memory
    "write_structures": False,
    "write_fasta": False,
    "checkpoint_path": mpnn_ckpt_dir
}
# Run sequence design on the RFD3-generated backbone
mpnn_model = MPNNInferenceEngine(**engine_config)

# Initialize RF3 inference engine
rf3_inference_engine = RF3InferenceEngine(ckpt_path=rf3_ckpt_dir, verbose=False)

print(f'!!! 记得，如果在每一步的out_dir 都赋值了，那么返回结果就会是None，不会保存在内存')
for design_length in range(length_min, length_max+1):
    print("="*30+f' ** 设计长度: {design_length} **'+"="*30+'\n')
    out_dir=f'{out_path}/{design_length}'
    # rfd3_out_dir = f'{out_dir}/rfd3'
    # mpnn_out_dir = f'{out_dir}/mpnn'
    # rf3_out_dir  = f'{out_dir}/rf3'

    # os.makedirs(rfd3_out_dir,exist_ok=True)
    # os.makedirs(mpnn_out_dir,exist_ok=True)
    # os.makedirs(rf3_out_dir,exist_ok=True)

    results_dir = out_dir
    # 生成时间戳用于文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # results_file = os.path.join(results_dir, f"pipeline_results_{timestamp}.json")
    results_file = os.path.join(results_dir, f"pipeline_results.json")

    # 用于存储所有结果的列表
    all_results = []
    print("\n"+"="*80+'\n')
    print(f'Section 1: All-Atom Generation with RFD3')
    # Configure RFD3 inference

    rfd3_seed = seed+design_length
    seed_everything(rfd3_seed)
    # Configure RFD3 inference(P23H-MD-1_E.pdb)
    config = RFD3InferenceConfig(
        specification={
            "dialect": 2,
            "input": input_pdb,
            # "contig": f"{design_length}-{design_length},/0,A1-348",  # P23H-MD-1_E.pdb
            "contig": f"{length_min}-{length_max},/0,A1-348",
            # "select_unfixed_sequence": "A20-35", # Converts selected indices in input to have unfixed sequence (inputs become atom14).
            # "length": design_length,
            },
        inference_sampler={
            "num_timesteps":200,
            "step_scale": 1.5,
        },
        diffusion_batch_size=diffusion_batch_size,  # Generate "diffusion_batch_size" structures per batch
        ckpt_path=rfd3_ckpt_dir,
        global_prefix=global_prefix,
        seed=rfd3_seed,
        dump_trajectories=True,
    )

    # Initialize engine and run generation
    model = RFD3InferenceEngine(**config)
    outputs = model.run(
        inputs=None,      # None for unconditional generation
        out_dir=None,     # None to return in memory (no file output)
        n_batches=rfd3_n_batches,      # Generate 1 batch
    )
    # 保存Section 1的结果
    # first_key = next(iter(outputs.keys()))
    # out  = outputs[first_key][0]
    # metrics=out.metadata['metrics']

    # 保存Section 1的结果
    section1_results = {}

    # Inspect RFD3 outputs and extract the generated structures
    for batch_idx, data in outputs.items():
        print(f"Batch {batch_idx}: {len(data)} structure(s)")
        for struct_idx, structure in enumerate(data):
            print(f"  Structure {struct_idx}: {type(structure).__name__}")
            print(f"    AtomArray shape: {structure.atom_array.array_length()}")
            
            # 保存Section 1的每个结果
            section1_key = f"batch_{batch_idx}_structure_{struct_idx}"
            section1_results[section1_key] = {
                "batch": batch_idx,
                "structure_index": struct_idx,
                "structure_type": type(structure).__name__,
                "atom_array_shape": structure.atom_array.array_length()
            }

    print("\n"+"="*80+'\n')
    print(f'Starting pipeline for all generated structures...')

    # 对Section 1生成的每个结构进行循环处理
    for batch_idx, data in outputs.items():
        batch_idx = str(batch_idx).replace("_","")
        for struct_idx, structure in enumerate(data):
            print(f"\nProcessing Batch {batch_idx}, Structure {struct_idx}")
            print("-" * 60)
            
            # 当前结构的唯一标识
            struct_id = f"batch_{batch_idx}_structure_{struct_idx}"
            
            # 初始化当前结构的结果字典
            struct_result = {
                "structure_id": struct_id,
                "section1": section1_results.get(struct_id, {}),
                "section2": {},
                "section3": {},
                "section4": {}
            }
            
            # 提取当前结构的原子阵列
            atom_array = structure.atom_array
            
            print(f'Section 2: Sequence Design with MPNN...')
            mpnn_seed =int(seed) + int(design_length) + int(batch_idx)+ int(struct_idx)

            # Configure per-input inference options
            # See mpnn.utils.inference.MPNN_PER_INPUT_INFERENCE_DEFAULTS for all options
            input_configs = [
                {
                    "batch_size": mpnn_batch_size if mpnn_batch_size is not None else 5, # Generate 10 sequences per structure
                    "remove_waters": True,
                    "number_of_batches": 1,
                    "seed": mpnn_seed,
                    "fixed_chains": ['B'],
                    "designed_chains": None
                }
            ]
            mpnn_outputs = mpnn_model.run(input_dicts=input_configs, atom_arrays=[atom_array])

            # 保存Section 2的结果
            designed_sequences = []
            
            print(f"Generated {len(mpnn_outputs)} designed sequences:\n")
            for i, item in enumerate(mpnn_outputs):
                res_starts = get_residue_starts(item.atom_array)
                mpnn_out = item.output_dict
                # Convert 3-letter codes to 1-letter using Biotite
                seq_1letter = ''.join(
                    ProteinSequence.convert_letter_3to1(res_name)
                    for res_name in item.atom_array.res_name[res_starts]
                )
                print(f"Sequence {i+1}: {seq_1letter}")
                
                # 保存序列信息
                designed_sequences.append({
                    "sequence_index": i,
                    "sequence": seq_1letter,
                    "atom_array_shape": item.atom_array.array_length(),
                    "output_dict": mpnn_out
                })
            
            struct_result["section2"] = {
                "total_designed_sequences": len(mpnn_outputs),
                "designed_sequences": designed_sequences,
                "first_sequence": designed_sequences[0]["sequence"] if designed_sequences else None
            }
            
            print("\n"+"="*80+'\n')
            print(f"Section 3: Structure Prediction with RF3")
            
            # 对MPNN设计的每个序列进行RF3预测（这里为了效率，可以选择只处理前几个，或全部处理）
            # 这里我们处理前3个序列作为示例
            rf3_results = []
            
            # 限制处理的序列数量以避免过长的运行时间
            max_sequences_to_process = min(10, len(mpnn_outputs))
            
            for seq_idx in range(max_sequences_to_process):
                print(f"\nProcessing designed sequence {seq_idx + 1}/{max_sequences_to_process}")
            
                # Create input from the MPNN-designed structure
                mpnn_atom_array = mpnn_outputs[seq_idx].atom_array
                example_id = f"{struct_id}_sequence_{seq_idx}"
                input_structure = InferenceInput.from_atom_array(mpnn_atom_array, example_id=example_id)
                rf3_outputs = rf3_inference_engine.run(inputs=input_structure, out_dir=None)

                # 提取预测结果
                rf3_sequence_results = []
                
                # 检查是否有输出
                if example_id in rf3_outputs and len(rf3_outputs[example_id]) > 0:
                    for model_idx, rf3_output in enumerate(rf3_outputs[example_id]):
                        # Summary confidences
                        summary = rf3_output.summary_confidences
                        
                        # Detailed per-atom/residue confidences
                        conf = rf3_output.confidences if rf3_output.confidences else {}
                        
                        # 保存RF3结果
                        rf3_result = {
                            "example_id": example_id,
                            "model_index": model_idx,
                            "atom_array_size": len(rf3_output.atom_array),
                            "atom_array": rf3_output.atom_array,  # 保存原子阵列以便后续使用
                            "summary_confidences": {
                                "overall_plddt": float(summary.get('overall_plddt', 0)),
                                "overall_pae": float(summary.get('overall_pae', 0)),
                                "overall_pde": float(summary.get('overall_pde', 0)),
                                "ptm": float(summary.get('ptm', 0)),
                                "iptm": float(summary.get('iptm', 0)) if summary.get('iptm') else None,
                                "ranking_score": float(summary.get('ranking_score', 0)),
                                "has_clash": bool(summary.get('has_clash', False))
                            },
                            "confidences": {
                                "atom_plddts_count": len(conf.get('atom_plddts', [])),
                                "first_10_atom_plddts": np.round(conf.get('atom_plddts', [])[:10], 2).tolist() if conf.get('atom_plddts') else [],
                                "pae_matrix_shape": f"{len(conf.get('pae', []))}x{len(conf.get('pae', [[]])[0])}" if conf.get('pae') else "0x0"
                            }
                        }
                        
                        rf3_sequence_results.append(rf3_result)
                else:
                    print(f"  Warning: No RF3 output for {example_id}")
                    rf3_sequence_results = []
                
                # 保存这个序列的RF3结果
                rf3_results.append({
                    "sequence_index": seq_idx,
                    "sequence": designed_sequences[seq_idx]["sequence"],
                    "rf3_models": rf3_sequence_results,
                    "best_model_summary": rf3_sequence_results[0]["summary_confidences"] if rf3_sequence_results else {}
                })
            
            struct_result["section3"] = {
                "processed_sequences": max_sequences_to_process,
                "total_sequences": len(mpnn_outputs),
                "rf3_predictions": rf3_results
            }
            
            print("\n"+"="*80+'\n')
            print(f'Section 4: Validation and Export for {struct_id}')
            
            # 对每个处理的序列进行验证
            validation_results = []
            
            for seq_idx in range(max_sequences_to_process):
                print(f"\nValidating designed sequence {seq_idx + 1}/{max_sequences_to_process}")
                
                # 获取RF3预测的最佳结构（第一个模型）
                if (seq_idx < len(rf3_results) and 
                    len(rf3_results[seq_idx]["rf3_models"]) > 0):
                    
                    # 从rf3_results中获取原子阵列
                    rf3_atom_array = rf3_results[seq_idx]["rf3_models"][0]["atom_array"]
                    
                    if rf3_atom_array is not None:
                        aa_refolded = rf3_atom_array
                        
                        # Get structures for comparison
                        aa_generated = atom_array  # Original RFD3 backbone
                        
                        # Filter to backbone atoms (N, CA, C, O)
                        bb_generated = aa_generated[np.isin(aa_generated.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
                        bb_refolded = aa_refolded[np.isin(aa_refolded.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
                        
                        # 确保两个结构有相同数量的原子用于比较
                        if len(bb_generated) == len(bb_refolded):
                            # Superimpose structures and calculate RMSD
                            bb_refolded_fitted, _ = superimpose(bb_generated, bb_refolded)
                            rmsd_value = rmsd(bb_generated, bb_refolded_fitted)
                            
                            # 解释RMSD值
                            if rmsd_value < 1.0:
                                interpretation = "Excellent"
                            elif rmsd_value < 2.0:
                                interpretation = "Good"
                            elif rmsd_value < 3.0:
                                interpretation = "Moderate"
                            else:
                                interpretation = "Poor"
                            
                            print(f"  Sequence {seq_idx + 1} Backbone RMSD: {rmsd_value:.2f} A ({interpretation})")
                            
                            # 保存验证结果
                            validation_results.append({
                                "sequence_index": seq_idx,
                                "sequence": designed_sequences[seq_idx]["sequence"],
                                "backbone_rmsd": float(rmsd_value),
                                "interpretation": interpretation,
                                "comparison_atoms": len(bb_generated)
                            })
                            
                            # Export structures to CIF format
                            struct_out_dir = os.path.join(out_dir, struct_id, f"sequence_{seq_idx}")
                            os.makedirs(struct_out_dir, exist_ok=True)
                            
                            to_cif_file(aa_generated, f"{struct_out_dir}/generated.cif")
                            to_cif_file(aa_refolded, f"{struct_out_dir}/refolded.cif")
                            
                            print(f"  Exported structures to: {struct_out_dir}")
                        else:
                            print(f"  Error: Backbone atom count mismatch ({len(bb_generated)} vs {len(bb_refolded)})")
                            validation_results.append({
                                "sequence_index": seq_idx,
                                "sequence": designed_sequences[seq_idx]["sequence"],
                                "backbone_rmsd": None,
                                "interpretation": "Error: Atom count mismatch",
                                "comparison_atoms": None
                            })
                    else:
                        print(f"  Error: No atom array available for sequence {seq_idx + 1}")
                        validation_results.append({
                            "sequence_index": seq_idx,
                            "sequence": designed_sequences[seq_idx]["sequence"],
                            "backbone_rmsd": None,
                            "interpretation": "Error: No RF3 prediction",
                            "comparison_atoms": None
                        })
                else:
                    print(f"  Warning: No RF3 models for sequence {seq_idx + 1}")
                    validation_results.append({
                        "sequence_index": seq_idx,
                        "sequence": designed_sequences[seq_idx]["sequence"],
                        "backbone_rmsd": None,
                        "interpretation": "Warning: No RF3 models",
                        "comparison_atoms": None
                    })
            
            # 计算平均RMSD（排除None值）
            valid_rmsds = [v["backbone_rmsd"] for v in validation_results if v["backbone_rmsd"] is not None]
            avg_rmsd = float(np.mean(valid_rmsds)) if valid_rmsds else 0.0
            
            struct_result["section4"] = {
                "validated_sequences": max_sequences_to_process,
                "validation_results": validation_results,
                "average_rmsd": avg_rmsd
            }
            
            # 将当前结构的结果添加到总结果列表
            all_results.append(struct_result)
            
            print(f"\nCompleted processing for {struct_id}")

    print("\n"+"="*80+'\n')
    print(f'Saving all results to JSON file: {results_file}')

    # 保存所有结果到JSON文件（注意：atom_array不能直接序列化为JSON，我们需要移除它）
    for struct_result in all_results:
        for rf3_pred in struct_result["section3"]["rf3_predictions"]:
            for model in rf3_pred["rf3_models"]:
                if "atom_array" in model:
                    del model["atom_array"]

    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    # 打印汇总信息
    print(f"\nPipeline Summary:")
    print(f"  Total structures processed: {len(all_results)}")
    print(f"  Total sequences designed: {sum([len(r['section2']['designed_sequences']) for r in all_results])}")
    print(f"  Total RF3 predictions: {sum([len(r['section3']['rf3_predictions']) for r in all_results])}")

    # 计算平均RMSD
    all_rmsds = []
    for struct_result in all_results:
        for validation in struct_result["section4"]["validation_results"]:
            if validation["backbone_rmsd"] is not None:
                all_rmsds.append(validation["backbone_rmsd"])

    if all_rmsds:
        avg_all_rmsd = np.mean(all_rmsds)
        print(f"  Average RMSD across all validated sequences: {avg_all_rmsd:.2f} Å")
    else:
        print(f"  No valid RMSD values to average")

    # 创建简化的汇总报告
    summary_report = {
        "pipeline_run_timestamp": timestamp,
        "input_pdb": input_pdb,
        "total_structures": len(all_results),
        "results_file": results_file,
        "structure_summaries": []
    }

    for struct_result in all_results:
        # 获取最佳RMSD
        valid_rmsds = [v["backbone_rmsd"] for v in struct_result["section4"]["validation_results"] 
                    if v["backbone_rmsd"] is not None]
        best_rmsd = min(valid_rmsds) if valid_rmsds else None
        
        struct_summary = {
            "structure_id": struct_result["structure_id"],
            "designed_sequences": struct_result["section2"]["total_designed_sequences"],
            "rf3_predictions": struct_result["section3"]["processed_sequences"],
            "best_rmsd": best_rmsd,
            "average_rmsd": struct_result["section4"]["average_rmsd"]
        }
        summary_report["structure_summaries"].append(struct_summary)

    # 保存汇总报告
    # summary_file = os.path.join(results_dir, f"pipeline_summary_{timestamp}.json")
    summary_file = os.path.join(results_dir, f"pipeline_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(summary_report, f, indent=2)

    print(f"\nSummary report saved to: {summary_file}")
    print("\n"+"="*80+'\n')
    print(f'Pipeline completed successfully!')

end = time.time()
elapsed_time = end - start  # 总耗时秒数
# 计算小时、分钟和剩余秒数
hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)
seconds = elapsed_time % 60  # 这里保留小数部分

print(f"代码执行耗时: {hours:02d}:{minutes:02d}:{seconds:06.3f}")
