import pandas as pd
import os
import csv
# ========================= 定义函数  =========================

def read_bed_pandas(file_path, ncols=5):
    """使用Pandas读取，适合大数据量"""
    columns = ['chrom', 'start', 'end', 'strand', 'id']
    df = pd.read_csv(file_path, sep='\t', header=None, comment='#', 
                     names=columns, dtype={'chrom': str})
    return df

def read_bed_pandas7(file_path, ncols,columns):
    """使用Pandas读取，适合大数据量"""
    # columns = ['chrom', 'start', 'end', 'strand', 'name',][:ncols]
    df = pd.read_csv(file_path, sep='\t', header=None, comment='#', 
                     names=columns, dtype={'chrom': str})
    return df

def read_tsv_pandas(file_path):
    """读取TSV到DataFrame，支持大文件分块读取"""
    # 常用参数：nrows=1000（读前N行）, chunksize=10000（分块读取）
    return pd.read_csv(file_path, sep='\t')

#定义计算和homolog及anchor相关的6个特征
def count_in_region(df, chrom_col, start_col, end_col, 
                            target_chrom, region_start, region_end):
    
    # 构建筛选条件
    condition = (
        (df[chrom_col] == target_chrom)&
        (df[start_col] >= region_start) &
        (df[end_col] <= region_end) 
    )
    # print(df[chrom_col])
    # &
    #     (df[start_col] >= region_start) &
    #     (df[end_col] <= region_end)
    # print(condition.sum())
    # 返回符合条件的行数
    return condition.sum()

def count_anchors_mapped(df,
                        chrom_col1,chrom_col2,start_col1, end_col1,  # 第一组坐标列
                         start_col2, end_col2,  # 第二组坐标列
                        target_chrom1, target_chrom2, region_start1, region_end1,  # 第一个目标区域
                         region_start2, region_end2): # 第二个目标区域  
                                       
    """
    统计同时位于两个指定区域的锚点数量

    参数：
    anchor_df : DataFrame - 包含双锚点信息的数据框
    chrom_col1/2 : str - 染色体列名（分别对应两个锚点）
    start_col1/2 : str - 起始位置列名
    end_col1/2 : str - 结束位置列名
    target_chrom1/2 : str - 目标染色体名称
    region_start1/2 : int - 目标区域起始位置
    region_end1/2 : int - 目标区域结束位置

    返回：
    int - 同时满足两个区域条件的行数
    """
    # 构建复合筛选条件
    condition = (
        (df[chrom_col1] == target_chrom1) &
        (df[start_col1] >= region_start1) &
        (df[end_col1] <= region_end1) &
        (df[chrom_col2] == target_chrom2) &
        (df[start_col2] >= region_start2) &
        (df[end_col2] <= region_end2)
    )
    # print(condition.sum())
    # print("%%%%%%%%%%%%%%%%%%%%%%%%%")
    # print(df[start_col1][0],"++",region_start1)
    # print(df[end_col1][0],"++",region_end1)
    # print(df[start_col2][0],"++",region_start2)
    # print(df[end_col2][0],"++",region_end2)

    return condition.sum()





#==============================gethomolog===============================
def read_mappings(mapping_file):
    """读取ID映射关系文件"""
    mappings = []
    with open(mapping_file) as f:
        next(f)  # 跳过标题行
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                mappings.append((parts[1], parts[2]))  # (nl_id, ag_id)
                # print(parts[1], parts[2])

    return mappings

def read_bed_file(bed_file):
    """读取BED格式文件，返回{ID: 前四列数据}的字典"""
    data = {}
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 5:
                data[cols[4]] = cols[:4]  # 使用第5列作为键
                # print(cols)
    return data

def gethomolog(mapping_file, bed1_file, bed2_file,columns):
    """主处理函数"""
    # 读取数据
    mappings = read_mappings(mapping_file)
    bed1_data = read_bed_file(bed1_file)
    bed2_data = read_bed_file(bed2_file)
    
    # 合并匹配结果
    result = []
    for nl_id, ag_id in mappings:
        if nl_id in bed1_data and ag_id in bed2_data:
            combined = bed1_data[nl_id] + bed2_data[ag_id]
            result.append(combined)
            # result.append('\t'.join(combined))
            # print(combined)
    # with open("homolog.bed", 'w') as f:
    #     for line in result:
    #         f.write(line + '\n')
    
    
    # # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(result,columns=columns)
    df['start1'] = df['start1'].astype(int)
    df['end1'] = df['end1'].astype(int)
    df['start2'] = df['start2'].astype(int)
    df['end2'] = df['end2'].astype(int)

    # # 保存为 CSV 文件（无索引）
    # df.to_csv("homolog.csv", index=False)

    return df



#========================= 相关性计算 =========================
def get_corr(nl_insect_anchor, nl_genome,insect_genome,nl_insect_pos):
    """
    计算数据框中指定列的相关性矩阵

    参数：
    df : DataFrame - 包含数据的数据框
    columns : list - 需要计算相关性的列名列表

    返回：
    DataFrame - 相关性矩阵
    """
    features_anchor = []
    nl_gene= nl_genome.iloc[0]
    insect_gene = insect_genome.iloc[0]  
    flanking_region = 1000000  # 1 Mbps

    #获取字典
    gene_dict = nl_genome.groupby(nl_genome.columns[4]).first().to_dict(orient='index')
    insect_gene_dict = insect_genome.groupby(insect_genome.columns[4]).first().to_dict(orient='index')
    # print(insect_gene_dict)

    for j in range(len(nl_insect_pos)):
    #    print("第",i,":",len(nl_insect_pos) ,"组")
        # if j==1:
        #     break  #qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq      
        nl_gene_name = nl_insect_pos.iloc[j, 0]


        # 查询时直接定位
        nl_gene = gene_dict.get(nl_gene_name)
        # for i in range(len(nl_genome)):
            
        #     # print(nl_genome.iloc[i, 4])
        #     if nl_gene_name == nl_genome.iloc[i, 4]:
        #         # print(nl_genome.iloc[i, 4])
        #         nl_gene= nl_genome.iloc[i]
        #         break
        
        insect_gene_name = nl_insect_pos.iloc[j, 1] 
        # 查询时直接定位
        insect_gene = insect_gene_dict.get(insect_gene_name)
        # for i in range(len(insect_genome)): 
        #     if insect_gene_name == insect_genome.iloc[i, 4]:
        #         insect_gene = insect_genome.iloc[i]     
        
        print("第",j,":",len(nl_insect_pos) ,"个") 
        # insect_gene = insect_genome.iloc[j]
        # print(insect_gene)
        # if not nl_gene.empty and not insect_gene.empty:
        if nl_gene is not None and insect_gene is not None:        
                nl_chrom = nl_gene['chrom']  # 使用方括号访问
                nl_start = nl_gene['start']
                nl_end = nl_gene['end']

                insect_chrom = insect_gene['chrom']
                insect_start = insect_gene['start']
                insect_end = insect_gene['end']
                m1u = count_in_region(nl_insect_anchor, 'chrom2', 'start2', 'end2', 
                                    nl_chrom, nl_start - flanking_region, nl_start)
                m1d = count_in_region(nl_insect_anchor, 'chrom2', 'start2', 'end2',
                                            nl_chrom, nl_end, nl_end + flanking_region)
                m1f = count_in_region(nl_insect_anchor, 'chrom2', 'start2', 'end2',
                                            nl_chrom, nl_start - flanking_region, nl_end + flanking_region)

                m2u = count_in_region(nl_insect_anchor, 'chrom1', 'start1', 'end1',
                                            insect_chrom, insect_start - flanking_region, insect_start)
                m2d = count_in_region(nl_insect_anchor, 'chrom1', 'start1', 'end1',
                                            insect_chrom, insect_end, insect_end + flanking_region)
                m2f = count_in_region(nl_insect_anchor, 'chrom1', 'start1', 'end1',
                                            insect_chrom, insect_start - flanking_region, insect_end + flanking_region)

                mu = count_anchors_mapped(nl_insect_anchor, 'chrom2', 'chrom1', 'start2', 'end2', 'start1', 'end1',
                                        nl_chrom, insect_chrom,
                                        nl_start - flanking_region, nl_start,
                                        insect_start - flanking_region, insect_start)
                md = count_anchors_mapped(nl_insect_anchor, 'chrom2', 'chrom1', 'start2', 'end2', 'start1', 'end1',
                                        nl_chrom, insect_chrom,
                                        nl_end, nl_end + flanking_region,
                                        insect_end, insect_end + flanking_region)
                mf = count_anchors_mapped(nl_insect_anchor, 'chrom2', 'chrom1', 'start2', 'end2', 'start1', 'end1',
                                        nl_chrom, insect_chrom,
                                        nl_start - flanking_region, nl_end + flanking_region,  # 修正原代码错误
                                        insect_start - flanking_region, insect_end + flanking_region)

                # 使用 Python 内置的 min 函数
                proportionscoreu = mu / min(m1u, m2u) if min(m1u, m2u) != 0 else 0
                proportionscored = md / min(m1d, m2d) if min(m1d, m2d) != 0 else 0
                proportionscoref = mf / min(m1f, m2f) if min(m1f, m2f) != 0 else 0

                features_anchor.append([
                    nl_gene_name,
                    insect_gene_name,  # 修正原代码中的 $$ 错误
                    mu, md, mf,
                    proportionscoreu, proportionscored, proportionscoref
                ])
                # print(nl_gene.name)
    return features_anchor


# ========================= 读取数据 =========================

kindss=["bm","cngd","ls","mifeng","sf","yachong"]
for kinds in kindss:
    root_dir = './训练集/nl_'+ kinds
    # 使用示例（读取BED6文件）
    nl_genome  = read_bed_pandas("./训练集/nlgenome.bed", ncols=5)
    path = os.path.join(root_dir, f'{kinds}genome.bed')
    insect_genome = read_bed_pandas(path,ncols=5)

    nl_lnc = read_bed_pandas("./预测集/nllncRNA.bed",ncols=5)

    path = os.path.join(root_dir, f'{kinds}lncRNA.bed')
    insect_lnc = read_bed_pandas(path,ncols=5)
    #读入锚点文件

    path = os.path.join(root_dir, f'nl_{kinds}_anchor.bed')
    nl_insect_anchor = read_bed_pandas7(path,ncols=7,
                                    columns = ['chrom1','start1','end1','chrom2','start2','end2','identity','extra1'])#1其他昆虫，2褐飞虱
    #获取homolog文件
    path = os.path.join(root_dir, f'nl_{kinds}_final2.tsv')
    nl_insect_pos = read_tsv_pandas(path)

    mapping_file = os.path.join(root_dir, f'nl_{kinds}_final.tsv')
    # mapping_file="./训练集/nl_bm/nl_bm_final.tsv"


    bed1_file="./训练集/nlgenome.bed"

    bed2_file= os.path.join(root_dir, f'{kinds}genome.bed')
    # bed2_file="./训练集/nl_bm/bmgenome.bed"

    nl_insect_homolog = gethomolog(mapping_file,bed1_file,bed2_file,columns = ['chrom2','start2','end2','strand2','chrom1','start1','end1','strand1'])  
    # nl_insect_anchor = read_bed_pandas7("./homolog.bed",ncols=7,
                                    #    columns = ['chrom2','start2','end2','strand2','chrom1','start1','end1','strand1'])#1其他昆虫，2褐飞虱                  

    # ======================== 计算阳性数据集和homolog及anchor相关的6个特征 ========================================
    # 计算和anchor相关的6个特征

    features_anchor = get_corr(nl_insect_anchor, nl_genome,insect_genome,nl_insect_pos)#nl_insect_anchor, nl_genome,insect_genome,nl_insect_pos
    # # flanking_region = 1000000  # 1 Mbps

    # print(features_anchor)

    # 写入输出文件
    import pandas as pd

    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(features_anchor, columns=[
        "nl_gene", "insect_gene",
        "mu", "md", "mf",
        "proportionscoreu", "proportionscored", "proportionscoref"
    ])

    # # 保存为 CSV 文件（无索引）
    savePath = os.path.join(root_dir, f'{kinds}train_features_anchor.csv')
    # df.to_csv("./训练集/nl_bm/train_features_anchor_pos.csv", index=False)
    df.to_csv(savePath, index=False)

    # 计算和homolog相关的6个特征==========================================================
    # print(nl_insect_homolog['start2'])

    # # 'chrom2', 'start2', 'end2'
    features_homolog = get_corr(nl_insect_homolog, nl_genome,insect_genome,nl_insect_pos)


    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(features_homolog, columns=[
        "nl_gene", "insect_gene",
        "mu", "md", "mf",
        "proportionscoreu", "proportionscored", "proportionscoref"
    ])

    # 保存为 CSV 文件（无索引）
    savePath = os.path.join(root_dir, f'{kinds}train_features_homolog_pos.csv')
    # df.to_csv("./训练集/nl_bm/train_features_homolog_pos.csv", index=False)
    df.to_csv(savePath, index=False)

    #数据整合
    nl_insect_pos_feature = [a + b for a, b in zip(features_anchor, features_homolog)]

    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(nl_insect_pos_feature, columns=[
        "nl_gene1", "insect_gene1",
        "mu1", "md1", "mf1",
        "proportionscoreu1", "proportionscored1", "proportionscoref1","nl_gene2", "insect_gene2",
        "mu2", "md2", "mf2",
        "proportionscoreu2", "proportionscored2", "proportionscoref2"
    ])

    # 保存为 CSV 文件（无索引）
    savePath = os.path.join(root_dir, f'{kinds}train_nl_insect_pos_feature.csv')
    # df.to_csv("./训练集/nl_bm/train_nl_insect_pos_feature.csv", index=False)
    df.to_csv(savePath, index=False)




    from itertools import product

    #============================== 计算阴性数据集和anchor及homolog相关的6个特征 ======================================================
    #代码逻辑同上，但输入文件不同
    #读取两物种阴性数据集的匹配id
    path = os.path.join(root_dir, f'nl_{kinds}_neg.csv')
    # nl_insect_neg = pd.read_csv("./训练集/nl_anwen_neg.csv")
    nl_insect_neg = pd.read_csv(path)
    # print(nl_insect_neg)

    features_anchor = get_corr(nl_insect_anchor, nl_genome,insect_genome,nl_insect_neg)#nl_insect_anchor, nl_genome,insect_genome,nl_insect_pos
    # # flanking_region = 1000000  # 1 Mbps

    # 写入输出文件
    import pandas as pd

    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(features_anchor, columns=[
        "nl_gene", "insect_gene",
        "mu", "md", "mf",
        "proportionscoreu", "proportionscored", "proportionscoref"
    ])

    # # 保存为 CSV 文件（无索引）
    path = os.path.join(root_dir, f'{kinds}train_features_anchor_neg.csv')
    # df.to_csv("train_features_anchor_neg.csv", index=False)
    df.to_csv(path, index=False)

    # 计算和homolog相关的6个特征==========================================================
    # print(nl_insect_homolog['start2'])

    # # 'chrom2', 'start2', 'end2'
    features_homolog = get_corr(nl_insect_homolog, nl_genome,insect_genome,nl_insect_neg)


    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(features_homolog, columns=[
        "nl_gene", "insect_gene",
        "mu", "md", "mf",
        "proportionscoreu", "proportionscored", "proportionscoref"
    ])

    # 保存为 CSV 文件（无索引）
    path = os.path.join(root_dir, f'{kinds}train_features_homolog_neg.csv')
    # df.to_csv("train_features_homolog_neg.csv", index=False)
    df.to_csv(path, index=False)
    #数据整合
    # nl_insect_neg_feature =
    nl_insect_neg_feature = [a + b for a, b in zip(features_anchor, features_homolog)]


    # 将列表转换为 DataFrame，并添加列名
    df = pd.DataFrame(nl_insect_neg_feature, columns=[
        "nl_gene1", "insect_gene1",
        "mu1", "md1", "mf1",
        "proportionscoreu1", "proportionscored1", "proportionscoref1","nl_gene2", "insect_gene2",
        "mu2", "md2", "mf2",
        "proportionscoreu2", "proportionscored2", "proportionscoref2"
    ])

    # 保存为 CSV 文件（无索引）
    path = os.path.join(root_dir, f'{kinds}train_nl_insect_neg_feature.csv')
    # df.to_csv("train_nl_insect_neg_feature.csv", index=False)
    df.to_csv(path, index=False)
