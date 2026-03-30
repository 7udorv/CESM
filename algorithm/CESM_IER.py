import copy
import time
import get_load
from ESPM_Col import CPM_Col
from ESPM_Col import BIESP_Col_incremental
import cProfile, pstats
import pandas as pd

def copes_miner_incre(count_dict, original_distances_dict, changed_distances_dict, r, min_prev, min_ei, min_overlap, original_files):
    #获取所有显著二元变化同位模式
    i = 0
    count_list = list(count_dict.keys())
    o_distances_list = list(original_distances_dict.keys())
    c_distances_list = list(changed_distances_dict.keys())
    ESPs = {}
    while True:
        print(i)
        if i == 0:
            instances1 = original_distances_dict[o_distances_list[i]]
            count1 = count_dict[count_list[i]]
            P1 = CPM_Col.cpm_col(instances1, r, min_prev, count1)
        instances2 = original_distances_dict[o_distances_list[i + 1]]
        count2 = count_dict[count_list[i + 1]]

        changed_instances = changed_distances_dict[c_distances_list[i]]

        Sig_ESP, P2 = BIESP_Col_incremental.Biesp_Col_incre(instances2, changed_instances, r, min_ei, min_prev, P1, count2)
        # print('end')
        ESPs[i] = Sig_ESP
        P1 = copy.deepcopy(P2)
        if i + 1 >= len(original_files) - 1:
            break
        else:
            i = i + 1

    if len(original_files) == 2:
        return

    n = 0
    P = {}
    start_time = time.perf_counter()
    while True:
        LS = set()
        if n == 0:
            C1 = [ESPs[n][cesp]['sequence'] for cesp in ESPs[n]]
        C2 = sorted([ESPs[n + 1][cesp]['sequence'] for cesp in ESPs[n + 1]])
        # 计算ISI
        # print('C1',C1)
        # print('C2',C2)
        for c1 in C1:
            flag_connectable = False
            c1_ep = c1.split('→')[-1]
            c1_ep_first = c1_ep.split(',')[0]
            for c2 in C2:
                c2_op = c2.split('→')[0]
                c2_op_first = c2_op.split(',')[0]
                if c1_ep_first < c2_op_first:
                    break
                else:
                    # print('c1_ep',c1_ep,'c2_op',c2_op)
                    if c1_ep == c2_op:
                        flag_connectable = True
                        flag_overlap = True
                        esp1_s = '→'.join(c1.rsplit('→', maxsplit=2)[-2:])
                        esp1 = [key for key in ESPs[n] if ESPs[n][key]['sequence'] == esp1_s][0]

                        esp2 = [key for key in ESPs[n + 1] if ESPs[n + 1][key]['sequence'] == c2][0]
                        # 首先找op的变化参与实例
                        op_dict = {}
                        ep_dict = {}
                        for f1 in c1_ep.split(','):
                            # 首先找esp1_s的实例
                            for f2 in esp1.split(','):
                                if f1 == f2.rstrip("+-"):
                                    ins_ep = ESPs[n][esp1][f2]
                                    if f1 == f2:  # 判断是否带变化特征
                                        ep_dict[f1] = ins_ep
                                    else:
                                        ins_ep = set(ins.rstrip("+-") for ins in ins_ep)
                                        ep_dict[f1] = ins_ep
                            # 其次找esp2_s的实例
                            for f3 in esp2.split(','):
                                if f1 == f3.rstrip("+-"):
                                    ins_op = ESPs[n + 1][esp2][f3]
                                    if f1 == f3:  # 判断是否带变化特征
                                        op_dict[f1] = ins_op
                                    else:
                                        ins_op = set(ins.rstrip("+-") for ins in ins_op)
                                        op_dict[f1] = ins_op

                            # 计算ISI
                            ISI = len(ep_dict[f1] & op_dict[f1]) / len(ep_dict[f1])
                            if ISI < min_overlap:
                                flag_overlap = False
                                break
                        if flag_overlap:
                            parts1 = c1.split("→")
                            parts2 = c2.split("→")
                            # 连接c1和c2
                            s = "→".join(parts1 + [parts2[-1]])
                            # print('s',s)
                            LS.add(s)
                        else:  # 不满足min_overlp
                            if n not in P:
                                P[n] = []
                            P[n].append(c1)
            if not flag_connectable:  # 如果不能连接
                if n not in P:
                    P[n] = []
                P[n].append(c1)

        C1 = C2 + list(LS)
        # 将C1排序
        C1 = sorted(C1)
        n = n + 1
        if n + 1 > len(original_files) - 2:
            P[n] = C1
            break
    # print(P)
    end_time = time.perf_counter()
    print(f"运行时间: {(end_time - start_time) * 1_000_000:.0f} 微秒")
    # print(P)
            #n+1表示终止时间

# r =800
# min_prev = 0.3
# min_ei = 0.3
# min_overlap = 0
#
# all_filepaths = ['../dataset/Beijing-converted/2017.csv','../dataset/Beijing-converted/2018.csv',
#                  '../dataset/Beijing-converted/2019.csv','../dataset/Beijing-converted/2020.csv']
# original_files = ['distances/origin_distances/Beijing_2018.csv', 'distances/origin_distances/Beijing_2019.csv',
#                   'distances/origin_distances/Beijing_2020.csv', 'distances/origin_distances/Beijing_2021.csv']
# changed_files = ['../distances/changed_distances/Beijing_2018-2019.csv', '../distances/changed_distances/Beijing_2019-2020.csv',
#                   '../distances/changed_distances/Beijing_2020-2021.csv']
#
#
# #获取所有count
# count_dict = {str(year): {} for year in range(2017, 2021)}
# for file_path in all_filepaths:
#     year = str(file_path.split('/')[-1].split('.')[0])  # 提取年份，例如 "2018"
#     count = get_load.get_count(file_path)
#     count_dict[year] = count # 存到字典里
#     # print(count)
#
#
# #获取所有距离列表
# original_distances_dict = {str(year): {} for year in range(2017, 2021)}
# for file_path in original_files:
#     year = file_path.split('_')[-1].split('.')[0]  # 提取年份，例如 "2018"
#     file_b = pd.read_csv(file_path)
#     instances = []
#     for index, row in file_b.iterrows():
#         instances.append([row['point1'], row['point2'], row['distance']])
#
#     original_distances_dict[year] = instances  # 存到字典里
#
# #获取所有变化距离
# changed_distances_dict = {f"{year}-{year+1}": {} for year in range(2017, 2020)}
# for file_path in changed_files:
#     year = file_path.split('_')[-1].split('.')[0]  # 提取年份，例如 "2018"
#     file_b = pd.read_csv(file_path)
#     instances = []
#     for index, row in file_b.iterrows():
#         instances.append([row['point1'], row['point2'], row['distance']])
#
#     changed_distances_dict[year] = instances  # 存到字典里
# # 把运行结果保存到临时文件
# cProfile.run(
#     'copes_miner_incre(count_dict, original_distances_dict, changed_distances_dict, r, min_prev, min_ei, min_overlap)',
#     filename='profile_data'
# )
#
# # 读取并显示前 10 个最耗时函数
# pstats.Stats('profile_data').sort_stats('tottime').print_stats(20)


