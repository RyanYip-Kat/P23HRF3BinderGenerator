import pandas as pd
import os
import shutil

def csv_to_excel(sheet_names, csv_files, output_excel):
    """
    将多个CSV文件写入同一个Excel文件的不同sheet中
    
    参数:
    sheet_names (list): 每个sheet的名称列表
    csv_files (list): 对应的CSV文件路径列表
    output_excel (str): 输出的Excel文件路径
    """
    outs =[]
    # 检查sheet名称和CSV文件数量是否一致
    if len(sheet_names) != len(csv_files):
        print("错误：sheet名称数量和CSV文件数量不匹配！")
        return
    
    # 定义需要删除的列名
    cols_to_drop = ['has_clash', 'backbone_rmsd', 'interpretation']
    # 创建Excel写入器，指定openpyxl引擎（支持多sheet）
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        for sheet_name, csv_file in zip(sheet_names, csv_files):
            try:
                # 检查CSV文件是否存在
                if not os.path.exists(csv_file):
                    raise FileNotFoundError(f"CSV文件不存在：{csv_file}")
                
                # 读取CSV文件
                df = pd.read_csv(csv_file)
                df = df.drop(columns=cols_to_drop, errors='ignore')
                print(f'{sheet_name} size : {df.shape}') 
                # 将DataFrame写入指定sheet
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"成功写入sheet: {sheet_name} (对应文件: {csv_file})")
                outs.append(df)
            
            except FileNotFoundError as e:
                print(f"警告：{e}，该sheet将被跳过")
            except Exception as e:
                print(f"警告：处理{csv_file}时出错 - {e}，该sheet将被跳过")
    return outs

# ------------------- 配置参数（根据你的实际情况修改） -------------------
if __name__ == "__main__":

    input_dir = "liuzhong/protein_binder/P23H_V2"
    combined_outdir= f"{input_dir}/final_data"
    os.makedirs(combined_outdir,exist_ok=True)
    sheet_names = ['batch1', 'batch2','batch3']

    csv_files = [f'{input_dir}/{sheet_name}/results_final.csv' for sheet_name in sheet_names]
    output_excel = f'{input_dir}/combined_results.xlsx'
    output_csv = f'{input_dir}/combined_results.csv'
    
    # 执行转换
    outs = csv_to_excel(sheet_names, csv_files, output_excel)
    print(f"\n转换完成！Excel文件已保存至：{output_excel}")
    df = pd.concat(outs,axis=0)
    # rfd3_file,rf3_file
    final_outs =[]
    for index,row in df.iterrows():
        rf3_file =  row['rf3_file']
        final_pdb_filename = os.path.basename(rf3_file)
        shutil.copyfile(rf3_file, combined_outdir+'/'+final_pdb_filename)
        row['final_pdb'] = combined_outdir+'/'+final_pdb_filename
        final_outs.append(row.to_dict())

    DF= pd.DataFrame(final_outs)
    DF = DF.drop(columns=["rfd3_file","rf3_file","rfd3_cif_file","rfd3_binder_seq","rf3_cif_file"], errors='ignore')
    DF.to_csv(output_csv,index=False)
    print(DF.shape)
