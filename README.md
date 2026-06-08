这是一份为您重新梳理和排版后的 `README.md` 文件。调整后的结构更加注重逻辑层次、可读性以及参数的清晰呈现，非常适合 GitHub/GitLab 等代码仓库的展示规范。

您可以直接复制以下 Markdown 内容使用：

```markdown
# RFdiffusion3 批量生成 Binder 自动化流程

本项目提供了一套完整的自动化流程，支持在给定目标蛋白质序列结构的基础上，基于 [RFdiffusion3](https://github.com/RosettaCommons/foundry) 批量生成指定长度范围的 Binder 结构，并自动完成各项评估指标的提取与分析。

完整流程划分为 5 个核心步骤，总执行顺序已封装在 `pipeline.sh` 脚本中供一键调用或分步执行。

---

## 环境配置

本流程的环境依赖与 RFdiffusion3 完全一致。请参考官方 GitHub 仓库完成环境搭建：
🔗 **[RosettaCommons/foundry](https://github.com/RosettaCommons/foundry)**

---

## 执行步骤

### 步骤 1：使用 RFdiffusion3 批量生成结构
调用核心大模型，根据设定的目标蛋白与参数，批量生成候选 Binder 结构。

```bash
python inference_loop.py

```

**配置参数说明（修改 `inference_loop.py` 中的变量）：**

* **控制参数**
* `seed`: 随机种子（如 `3872541`），用于确保结果的计算可复现性。
* `length_min` / `length_max`: 生成 Binder 的最小与最大长度范围（例如 `70` 与 `86`）。
* `diffusion_batch_size`: 每次 RFD3 扩散生成的结构数量。
* `rfd3_n_batches`: 运行 RFD3 的总批次数。
* `mpnn_batch_size`: 提供给 ProteinMPNN 处理的序列数量。


* **路径参数**
* `input_pdb`: 目标蛋白的 PDB 文件路径（如 `"liuzhong/P23H_opsin_chainA.pdb"`）。
* `global_prefix`: 生成结构文件的前缀标识（如 `'p23h'`）。
* `out_path`: 当前批次结果的保存目录（如 `"/your/path/batch1"`）。



> **💡 结构产出说明：**
> * 总生成文件数量 = `diffusion_batch_size` × `rfd3_n_batches` × `mpnn_batch_size` × (`length_max` - `length_min`)
> * **进阶注意（Contig 设置）：** 在 `RFD3InferenceConfig` 中，`contig` 参数格式需设为 `"{length_min}-{length_max},/0,A1-348"`。其中 `A` 代表 `input_pdb` 中序列所在的链标识，`348` 为目标链的序列总长度。
> 
> 

---

### 步骤 2：单批次结果处理与指标提取

对步骤 1 生成的独立批次（如 batch1, batch2 等）进行文件解析与基础指标统计。

```bash
python process_design_results.py 

```

**配置参数说明：**

* `batch_id`: 当前处理的批次标识（如 `"batch1"`）。
* `original_pdb_file`: 原始目标 PDB 路径（需与步骤 1 的 `input_pdb` 保持一致）。
* `result_dir`: 当前批次的读取路径（对应步骤 1 的 `out_path`）。
* `result_csv_file`: 单批次结果指标的输出路径（如 `f"{result_dir}/results.csv"`）。
* `seq_ranges`: Binder 生成的长度区间列表（如 `[70, 86]`，需与步骤 1 保持一致）。

---

### 步骤 3：多批次结果统一汇总

将各个独立 Batch 的解析结果与 PDB 结构进行整合，便于后续环节的批量打分。

```bash
python combined_final_results.py

```

**配置参数说明：**

* `input_dir`: 包含所有 Batch 文件夹的根级目录（同步骤 1 的上级路径）。
* `sheet_names`: 需要汇总的独立批次文件夹名称列表（如 `['batch1', 'batch2', 'batch3']`）。
* `combined_outdir`: 汇总后统一存放结构与数据的输出目录（如 `f"{input_dir}/final_data"`）。

---

### 步骤 4：Rosetta InterfaceAnalyzer 打分

调用 Rosetta 工具对汇总后的 Binder 结构进行界面结合能（Interface Energy）等物理化学指标计算。

```bash
bash run_rosetta_InterfaceAnalyzer.sh /path/to/final_data

```

**参数说明：**

* `/path/to/final_data`: 必须指向步骤 3 中设定好的 `combined_outdir` 汇总目录。

---

### 步骤 5：最终结果合并

将步骤 3 提取的基础结构特征与步骤 4 计算的 Rosetta 结合能指标进行最终对齐与合并，输出完整数据表。

```bash
python merge_rosetta_result.py

```

**配置参数说明：**

* `step3_csv_file`: 步骤 3 产出的结构基础指标文件（如 `"liuzhong/protein_binder/P23H_V2/combined_results.csv"`）。
* `step4_csv_file`: 步骤 4 产出的 Rosetta 打分结果文件（如 `"liuzhong/protein_binder/P23H_V2/final_data/InterfaceAnalyzer_output.csv"`）。

