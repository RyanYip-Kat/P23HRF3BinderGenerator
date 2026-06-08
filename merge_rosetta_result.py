import pandas as pd
import re

def merge_dfs_by_pdb_description(df1, df2, pdb_col='final_pdb', desc_col='description', merge_how='outer'):
    """
    根据df1的final_pdb列和df2的description列合并两个DataFrame
    参数说明：
        df1: 包含final_pdb列的DataFrame（原表格1）
        df2: 包含description列的DataFrame（原表格2）
        pdb_col: df1中pdb路径的列名，默认'final_pdb'
        desc_col: df2中描述字段的列名，默认'description'
        merge_how: 合并方式，可选'inner/left/right/outer'，默认'outer'（保留所有行）
    返回值：
        合并后的DataFrame
    注意：
        1. 函数会复制原数据，不会修改输入的df1/df2
        2. 空值会被标记为None，不参与匹配
    """
    # 复制数据，避免修改原表格
    df1_copy = df1.copy()
    df2_copy = df2.copy()

    # -------------------------- 处理df1的final_pdb列 --------------------------
    def extract_core_from_pdb(pdb_path):
        """从pdb路径中提取核心标识"""
        if pd.isna(pdb_path):  # 处理空值
            return None
        # 拆分路径取最后一个文件名（如从路径中提取batch1_xxx.pdb）
        filename = pdb_path.split('/')[-1]
        # 去掉.pdb后缀，得到核心标识
        core_id = filename.replace('.pdb', '')
        return core_id

    # 新增临时列存储核心标识
    df1_copy['_temp_core_id'] = df1_copy[pdb_col].apply(extract_core_from_pdb)

    # -------------------------- 处理df2的description列 --------------------------
    def extract_core_from_desc(desc):
        """从description中提取核心标识（去掉末尾_数字后缀）"""
        if pd.isna(desc):  # 处理空值
            return None
        # 正则匹配：去掉末尾的"_+数字"（如_0001、_123）
        core_id = re.sub(r'_\d+$', '', desc)
        return core_id

    # 新增临时列存储核心标识
    df2_copy['_temp_core_id'] = df2_copy[desc_col].apply(extract_core_from_desc)

    # -------------------------- 合并两个表格 --------------------------
    merged_df = pd.merge(
        df1_copy, df2_copy,
        on='_temp_core_id',       # 基于核心标识合并
        how=merge_how,            # 合并方式
        suffixes=('_df1', '_df2') # 同名列添加后缀区分
    )

    # 删除临时核心标识列（如需保留可注释此行）
    merged_df = merged_df.drop('_temp_core_id', axis=1)

    return merged_df

# -------------------------- 测试用例（可直接运行） --------------------------
if __name__ == "__main__":

    step3_csv_file = "liuzhong/protein_binder/P23H_V2/combined_results.csv"
    step4_csv_file = "liuzhong/protein_binder/P23H_V2/final_data/InterfaceAnalyzer_output.csv"
    df1=pd.read_csv(step3_csv_file)
    df2=pd.read_csv(step4_csv_file)

    df1_use_col = 'final_pdb'
    df2_use_col = 'description'
    # 调用函数合并
    merged_result = merge_dfs_by_pdb_description(df1, df2,df1_use_col,df2_use_col,merge_how='inner')   
    # 打印结果
    merged_result.to_csv("liuzhong/protein_binder/P23H_V2/final_result_data.csv",index=False)
    print(merged_result.shape)
