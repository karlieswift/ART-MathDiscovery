import pandas as pd
import os
import sys
import argparse
import time
import numpy as np

standard_aa = list("ACGT")
#standard_aa = list("ACDEFGHIKLMNPQRSTVWY")  # 标准20种氨基酸
VALID_CHARS = "ACGT"

def calculate_values(sequence):
    nucleotides = 'ACGT'
    n = {nucl: 0 for nucl in nucleotides}
    mu = {nucl: 0 for nucl in nucleotides}
    D2 = {nucl: 0 for nucl in nucleotides}

    for nucl in nucleotides:
        n[nucl] = sequence.count(nucl)

    for nucl in nucleotides:
        positions = [i + 1 for i, char in enumerate(sequence) if char == nucl]
        if positions:
            mu[nucl] = sum(positions) / n[nucl]

    total_length = len(sequence)
    for nucl in nucleotides:
        positions = [i + 1 for i, char in enumerate(sequence) if char == nucl]
        if positions:
            D2[nucl] = sum((pos - mu[nucl]) ** 2 for pos in positions) / total_length / n[nucl]

    return n, mu, D2

# k-mer
# return dict[tag] = value
def calculate_all_NV_Kmer_sliding_window(sequence, mu, n, k=2):
    tmp_key_dict = dict()
    time_start = time.perf_counter()

    if k > len(sequence):
        return dict()
    # init
    covariances = dict()
    seq_mu_list = np.array([mu[aa] for aa in sequence], dtype=object)

    # sigma(index:1:n): (index - mu1)(index - mu2)(index - mu3)

    # enumerate
    for i in range(len(sequence)-k+1):

        tag = sequence[i:i+k]

        # build tmp record set of calc-ed index of sth tag
        if tag not in tmp_key_dict:
            tmp_key_dict[tag] = set()

        # calc tags' all feature value, and dont sum added value
        for j in range(i, i+k):
            if j not in tmp_key_dict[tag]:
                # record
                tmp_key_dict[tag].add(j)
                # calc
                #res = (j + 1 - seq_mu_list[i:i + k]).prod()
                res = abs((j + 1 - seq_mu_list[i:i + k])).prod()
                # add
                if tag in covariances:
                    covariances[tag] += res
                else:
                    covariances[tag] = res

    time_tag_enum = time.perf_counter()
    #print("new slide no sep")
    #print(covariances)
    for key_str in covariances:
        # [a, b, c, d, ...]
        under_n = np.array([n[a] for a in key_str], dtype=object).prod()
        covariances[key_str] = covariances[key_str] / under_n

    time_update = time.perf_counter()
    #print("[new slide]", "===" * 9)
    #print([(k, covariances[k]) for k in sorted(covariances)])

    all_time = time_update - time_start
    #print(">>Total time: %.4fs" % (time_update - time_start))

    return covariances

def calc_specific_NV_Kmer_sliding_window(sequence, mu, n, input_tag):
    """计算特定k-mer的滑动窗口特征"""
    tmp_key_dict = dict()
    time_start = time.perf_counter()

    # 初始化
    covariances = {tag: 0 for tag in input_tag}
    seq_mu_list = np.array([mu[aa] for aa in sequence], dtype=object)

    # 枚举所有标签
    for tag in input_tag:
        if tag not in tmp_key_dict:
            tmp_key_dict[tag] = set()

        start_idx = 0
        while start_idx < len(sequence):
            if tag not in sequence[start_idx:]:
                break

            find_idx = sequence[start_idx:].index(tag)
            start_idx = find_idx + start_idx  # 跳转到找到的位置

            for j in range(start_idx, start_idx + len(tag)):
                if j not in tmp_key_dict[tag]:
                    # 记录
                    tmp_key_dict[tag].add(j)
                    # 计算
                    #res = (j + 1 - seq_mu_list[start_idx: start_idx + len(tag)]).prod()
                    res = abs((j + 1 - seq_mu_list[start_idx: start_idx + len(tag)])).prod() #abs better???
                    # 累加
                    covariances[tag] += res

            start_idx += 1

    time_tag_enum = time.perf_counter()

    # 标准化 - 修复除以零的问题
    for key_str in covariances:
        try:
            # 检查k-mer中的每个氨基酸是否都在n字典中
            under_n_components = []
            for a in key_str:
                if a in n and n[a] != 0:
                    under_n_components.append(n[a])
                else:
                    # 如果某个氨基酸不存在或计数为0，设置under_n为0,计算covariances检查
                    under_n_components.append(0)
                    #print(f"Warning: Amino acid '{a}' in k-mer '{key_str}' not found in sequence or count is zero")
            
            under_n = np.array(under_n_components, dtype=object).prod()
            
            if under_n != 0:
                covariances[key_str] = covariances[key_str] / under_n
            else:
                covariances[key_str] = 0
                #print(f"Warning: under_n is zero for k-mer '{key_str}', setting covariance to 0")
            if covariances[key_str] > 1E+20 :
                print(f"Warning: covariances is too big for k-mer '{key_str}'")

        except Exception as e:
            print(f"Error processing k-mer '{key_str}': {e}")
            covariances[key_str] = 0

    time_update = time.perf_counter()
    #print(f"[new slide specific] Total time: {time_update - time_start:.4f}s")

    return covariances

def load_vocab(vocab_file):
    """加载词汇表文件"""
    vocab = []
    try:
        with open(vocab_file, 'r') as f:
            next(f)  # 跳过标题行
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    token = parts[0].strip()
                    if token:  # 确保token非空
                        vocab.append(token)
        print(f"Loaded {len(vocab)} tokens from {vocab_file}")
        return vocab
    except Exception as e:
        print(f"Error loading vocab file: {e}")
        return []

def process_seq_(sequence, vocab_tokens=None):
    """处理序列，计算NV特征"""
    n, mu, D2 = calculate_values(sequence)
    NV_n_mu = list(n.values()) + list(mu.values())  # 4 + 4 = 8个特征

    # 如果提供了词汇表，计算额外的k-mer特征
    if vocab_tokens:
        kmer_features = calc_specific_NV_Kmer_sliding_window(sequence, mu, n, vocab_tokens)
        # 按词汇表顺序排列特征值
        kmer_values = [kmer_features.get(token, 0) for token in vocab_tokens]
        NV_final = NV_n_mu + kmer_values
        #return kmer_values
        return NV_final
    else:
        # 如果没有词汇表，返回空列表或其他默认值
        return []

def generate_kmers(alphabet, k):
    """生成指定长度的所有kmer"""
    if k == 1:
        return alphabet
    else:
        kmers = []
        for mer in generate_kmers(alphabet, k-1):
            for aa in alphabet:
                kmers.append(mer + aa)
        return kmers

def process_seq(sequence, vocab_tokens=None, K=6):
    """处理序列，计算NV特征"""
    n, mu, D2 = calculate_values(sequence)
    NV_n_mu = list(n.values()) + list(mu.values())  # 4 + 4 = 8个特征

    # 生成所有kmer并计算特征
    all_kmers_list = []
    NV_kmer = []
    '''    
    for k in range(2, K+1):  # 生成2,3,4,5-mer
        # 为每个k值单独计算kmer特征
        NV_dict = calculate_all_NV_Kmer_sliding_window(sequence, mu, n, k)
        kmers = generate_kmers(standard_aa, k)
        all_kmers_list.extend(kmers)
        
        # 按固定顺序获取当前k值的kmer特征值
        current_kmer_features = [NV_dict.get(kmer, 0) for kmer in kmers]
        NV_kmer.extend(current_kmer_features)

    # 基础特征 + kmer特征
    NV = NV_n_mu + NV_kmer

    # 计算预期特征数量
    n_mu_features = len(standard_aa) * 2  # n和mu各4个
    kmer_features = sum(len(generate_kmers(standard_aa, k)) for k in range(2, K+1))
    total_expected = n_mu_features + kmer_features

    #print(f"最终特征长度: {len(NV)} (预期: {total_expected})")
    #print(f"基础特征(n+mu): {n_mu_features}, kmer特征: {kmer_features}")

    # 验证特征数量
    if len(NV) != total_expected:
        print(f"警告: 特征数量不匹配! 实际: {len(NV)}, 预期: {total_expected}")

    # 如果提供了词汇表，计算额外的k-mer特征
    '''
    if vocab_tokens:
        kmer_features = calc_specific_NV_Kmer_sliding_window(sequence, mu, n, vocab_tokens)
        # 按词汇表顺序排列特征值
        kmer_values = [kmer_features.get(token, 0) for token in vocab_tokens]
        #NV_final = NV_kmer + kmer_values
        NV_final = NV_n_mu + kmer_values
        return NV_final

    #return NV

def clean_sequence(index, sequence, valid_chars):
    """Remove invalid characters from the sequence and print them out."""
    # 1. 转换为大写并替换U为T
    sequence = sequence.upper().replace('U', 'T')
    
    cleaned_sequence = ''
    removed_chars = []
    
    for char in sequence:
        if char not in valid_chars:
            removed_chars.append(char)
        else:
            cleaned_sequence += char

    #if removed_chars:
    #    print(f"Invalid characters found at {index} and removed: {set(removed_chars)}")

    return cleaned_sequence

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='Process protein sequences and generate NV features')
    parser.add_argument('--input', required=True, help='Input CSV file path')
    parser.add_argument('--output', required=True, help='Output CSV file path')
    parser.add_argument('--vocab', help='Vocabulary file path for k-mer features')
    args = parser.parse_args()

    # 加载词汇表
    vocab_tokens = None
    if args.vocab:
        vocab_tokens = load_vocab(args.vocab)
        if not vocab_tokens:
            print("Warning: No vocabulary tokens loaded, skipping k-mer features")
        else:
            print(f"First 10 tokens: {vocab_tokens[:10]}")
            print(f"Last 10 tokens: {vocab_tokens[-10:]}")

    # 读取数据
    print(f"Loading data from {args.input}")
    reshaped_data = pd.read_csv(args.input)
    print(reshaped_data.head())
    print(f"Data shape: {reshaped_data.shape}")
    print(f"Valid characters: {VALID_CHARS}")

    # 处理每一行，生成NV特征
    nv_data = []
    total_rows = len(reshaped_data)

    for index, row in reshaped_data.iterrows():
        if index % 10000 == 0:
            print(f"Processing row {index}/{total_rows}")
        
        # 构建序列 - 根据你的数据结构调整这里
        sequence = row['Sequence']
        #sequence = row['sequence']
        cleaned_sequence = clean_sequence(index, str(sequence), VALID_CHARS)
        
        # 检查清理后的序列长度
        if len(cleaned_sequence) == 0:
            print(f"Warning: All characters removed from sequence at row {index}, skipping")
            continue
            
        try:
            NV = process_seq(cleaned_sequence, vocab_tokens)
            nv_data.append(NV)
        except Exception as e:
            print(f"Error processing sequence at row {index}: {e}")

    # 自动生成列名 - 简化版本
    if nv_data:
        actual_feature_count = len(nv_data[0])
        nv_columns = ['V' + str(i + 1) for i in range(actual_feature_count)]
        print(f"自动生成 {actual_feature_count} 个特征列")
    
        nv_df = pd.DataFrame(nv_data, columns=nv_columns)
    else:
        print("错误: 没有特征数据")
        nv_df = pd.DataFrame()
    
    print(nv_df.head())
    print(f"NV features shape: {nv_df.shape}")

    # 合并原始数据和NV特征
    processed_data = pd.concat([reshaped_data.reset_index(drop=True), nv_df], axis=1)
    print(processed_data.head())
    print(f"Processed data shape: {processed_data.shape}")

    # 保存结果
    print(f"Saving processed data to {args.output}")
    processed_data.to_csv(args.output, index=False)
    print("Processing completed successfully")

if __name__ == "__main__":
    main()
