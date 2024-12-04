import numpy as np
import pandas as pd
from twobitreader import TwoBitFile

genome = TwoBitFile("hg38.2bit")

def get_sequence(array):

    chrom_from_array = array[0]
    start = array[1]
    end = array[2]
    
    # 处理可能的异常输入
    try:
        chrom = genome[chrom_from_array]
        seq = chrom[start:end]
    except KeyError:
        raise ValueError(f"Chromosome {chrom_from_array} not found in genome")
    return seq

def preprocess(filename):

    name = filename.split('.')[0]
    
    # 读取 CSV 文件并选择所需列
    df = pd.read_csv(filename)
    
    # 选取前三列（染色体，起始位置，结束位置）
    df_selected = df.iloc[:, [0, 1, 2]]
    
    # 转换为 NumPy 数组以便应用向量化的操作
    df_selected_numpy = df_selected.to_numpy()
    
    # 使用 apply_along_axis 获取序列并将其转换为大写
    df_selected_numpy_seq = np.apply_along_axis(get_sequence, 1, df_selected_numpy)
    df_selected_numpy_seq = np.char.upper(df_selected_numpy_seq)
    
    # 将获取到的序列加入到原数据框
    df['seq'] = df_selected_numpy_seq
    
    # 选择并重命名所需的列
    dataset = df[['seq', 'peak_k562:ctcf']].rename(columns={'peak_k562:ctcf': 'result'})
    
    # 保存结果到新文件中
    dataset.to_csv(f"{name}_seq_result.csv", index=False)

if __name__ == "__main__":
    # 预处理训练和测试集
    preprocess("K562_CTCF_train.csv")
    preprocess("K562_CTCF_test.csv")


    
    