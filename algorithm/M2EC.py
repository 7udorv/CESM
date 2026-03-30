from copy import deepcopy
from itertools import combinations
import copy
import pandas as pd
from ESPM_Col import CPM_Col
import get_load
from ESPM_Col import CPM_Col

count_col = 0

def Biesp_Col(instances2, changed_instances, r, min_ei, min_prev, P1, count2):

    count_candidates = 0
    P2 = CPM_Col.cpm_col(instances2, r, min_prev, count2)
    # 获取所有变化点及其邻居的邻居
    # print('P1', P1)
    # print('P2',P2)
    cgroupN, neighbors = CPM_Col.get_groupN(changed_instances, r)
    # 根据neighbors生成所有二阶变化同位模式,其中包括不包含变化特征的模式
    esp2 = gen_esp2(neighbors)
    # print(len(esp2))
    k = 2
    All_Eobj = {}
    All_Eobj[k] = esp2
    next_candidates = esp2.keys()
    Sig_Esp = {}
    while True:
        cesps = gen_cesp(list(next_candidates), k)
        count_candidates += len(cesps)
        # print('cesps', cesps)
        next_candidates = copy.deepcopy(cesps)
        ESPs_now = {}
        for cesp in cesps:
            #首先检查cesp中是否包含变化特征，不包含则直接加入候选
            if '+' in cesp or '-' in cesp: #包含变化特征
                #检查原模式是否频繁，对于在P中的，直接认为op频繁，对于不在P中的，需要进一步检验频繁性
                flag, op, ep = is_OpPrevalent(cesp, P1)
                #记录已找到的变化op
                if flag:
                    # 根据k-1阶现有的演化参与实例剪枝
                    # for key in All_Eobj[k].keys():
                    #     flag_suer = True
                    #     if set(key.split(',')).issubset(set(cesp.split(','))):   #如果是子集
                    #         #前面已经判断过op
                    #         for f in key.split(','):
                    #             if f.rstrip('+-') in P1[len(op.split(','))][op]:
                    #                 super_UER = len(All_Eobj[k][key][f])/len(P1[len(op.split(','))][op][f.rstrip('+-')])
                    #                 if super_UER < min_ei:
                    #                     # print("super_pruned")
                    #                     flag_suer = False
                    #                     break
                    #         if not flag_suer:
                    #             break    #直接跳出第一层循环
                    # if not flag_suer:    #cesp的super_UER不满足条件
                    #     len_ep = len(ep.split(','))
                    #     if len_ep not in P2 or ep not in P2[len_ep]:   #如果ep不频繁
                    #         next_candidates.remove(cesp)
                    #     continue   #跳过

                    #获取cesp的所有演化参与实例
                    eobj = find_eobj(cesp, k, All_Eobj, P1, op, min_ei, cgroupN)
                    # print('eobj', eobj)
                    if eobj == {}: #表示没有参与实例，删除候选
                        next_candidates.remove(cesp)
                    if len(eobj) != 0:  #eobj不为空
                        ESPs_now[cesp] = eobj
                        #检查cesp是否满足显著阈值
                        for f in cesp.split(','):
                            flag_ei = True
                            if f.rstrip('+-') in op.split(','):
                                ER = len(eobj[f])/len(P1[len(op.split(','))][op][f.rstrip('+-')])
                                # op_dict[f.rstrip('+-')] = eobj[f]
                                if ER < min_ei:
                                    flag_ei = False
                                    break
                            # if f.rstrip('+-') in ep.split(','):
                            #     ep_dict[f.rstrip('+-')] = eobj[f]
                        if flag_ei:  #保存所有op频繁，并且满足min_ei的cesp
                            if cesp not in Sig_Esp:
                                Sig_Esp[cesp] = eobj
                                Sig_Esp[cesp]['sequence'] = f'{op}→{ep}'
                                # Sig_Esp[cesp]['op_dict'] = op_dict
                                # Sig_Esp[cesp]['ep_dict'] = ep_dict

                else:#如果op不频繁，从候选中删除cesp
                    next_candidates.remove(cesp)

            else:
                # 如果不包含变化特征,直接检查cesp是否频繁，如果不频繁，连接后的cesp的op也必定不频繁
                if len(cesp.split(',')) in P1 and cesp in P1[len(cesp.split(','))]:
                    flag = True
                else:
                    flag = False
                    next_candidates.remove(cesp)
                if flag:
                    #获取cesp的所有演化参与实例
                    eobj = find_eobj(cesp, k, All_Eobj, P1, cesp, min_ei, cgroupN)
                    if eobj:  #eobj不为空
                        ESPs_now[cesp] = eobj
        # print('S',Sig_Esp)
        # 根据Sig_Esp中所有模式获取P2
        for cesp in list(next_candidates):
            cesp_list = cesp.split(',')
            op, ep = [], []
            # 首先寻找op和ep
            for f in cesp_list:
                if f[-1] == '+':
                    ep.append(f.rstrip('+'))
                elif f[-1] == '-':
                    op.append(f.rstrip('-'))
                else:
                    op.append(f)
                    ep.append(f)
            ep = ','.join(ep)
            # 在P2中则频繁
            if len(ep.split(',')) not in P2 or ep not in P2[len(ep.split(','))]:
                if cesp in Sig_Esp.keys():
                    del Sig_Esp[cesp]
                next_candidates.remove(cesp)

        if len(next_candidates) == 0:
            break

        k = k + 1
        All_Eobj[k] = ESPs_now

        if len(next_candidates) == 0:
            break
    print(count_col)
    print('count_candidates',count_candidates)
    return Sig_Esp, P2



def gen_esp2(neighbors):
    candidates = {}
    for neighbor in neighbors:
        pt1 = neighbor[0]
        f1  = pt1.split('.')[0]
        pt2 = neighbor[1]
        f2 = pt2.split('.')[0]
        if pt1[-1] == '+' or pt1[-1] == '-':
            f1 = f1 + pt1[-1]
        if pt2[-1] == '+' or pt2[-1] == '-':
            f2 = f2 + pt2[-1]
        features = f1 +','+ f2
        if features not in candidates.keys():
            candidates[features] = {}
            candidates[features][f1] = set()
            candidates[features][f2] = set()
        candidates[features][f1].add(pt1)
        candidates[features][f2].add(pt2)
    return candidates


def gen_cesp(esp_list, k):
    esp_list = sorted(esp_list, key = lambda x: x[0])
    cesp = set()
    for i in range(0, len(esp_list)-1):
        for j in range(i+1, len(esp_list)):
            esp1 = esp_list[i].split(',')
            esp2 = esp_list[j].split(',')
            if esp1[0].rstrip('+-') <esp2[0].rstrip('+-'):
                break
            elif esp1[0:k-1] == esp2[0:k-1]:
                if esp1[-1].rstrip('+-') != esp2[-1].rstrip('+-'):  #要求最后一个元素的特征不能相同
                    if esp1[-1] < esp2[-1]:
                        c = esp1 + [esp2[-1]]
                    else:
                        c = esp2 + [esp1[-1]]
                    c = ','.join(c)
                    cesp.add(c)
    return cesp

def is_OpPrevalent(cesp, P):
    cesp_list = cesp.split(',')
    op, ep = [], []
    #首先寻找op和ep
    for f in cesp_list:
        if f[-1] == '+':
            ep.append(f.rstrip('+'))
        elif f[-1] == '-':
            op.append(f.rstrip('-'))
        else:
            op.append(f)
            ep.append(f)
    len_op = len(op)
    op  = ','.join(op)
    ep = ','.join(ep)
    #判断op的频繁性
    flag = False
    #若len_op为1，跳过
    if len_op <= 1:
        return flag, op, ep
    else:
        if len(op.split(',')) in P and op in P[len_op]:
            flag = True
        else:
            flag = False
    return flag, op, ep

def find_eobj(cesp, k, All_Eobj, P, op, min_ei, cgroupN):
    global count_col
    #获取cesp的候选演化参与实例
    # 获取候选参与实例
    co_eobj = {}
    eobj = {}   #用于存放所有演化模式的演化参与实例
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = set()
                co_eobj[f] = co_eobj[f] | All_Eobj[k][sub][f]



    #获取op的obj
    len_op = len(op.split(','))

    if len(co_eobj) != len(cesp.split(',')):  #说明缺失的特征没有演化参与实例
        return {}
    #     #从k-1阶找到缺失的f的演化参与实例
    #     print('All_Eobj',All_Eobj)
    #     print('co_eobj', co_eobj)
    #     print('cesp', cesp)
    #     missing_features = set(cesp.split(',')) - set(co_eobj.keys())
    #     print('missing_features', missing_features)

    # # 计算演化参与率上界(只有演化同位模式需要)
    # if '+' in cesp or '-' in cesp:
    #     for f in cesp.split(','):
    #         if f.rstrip('+-') in op.split(','):
    #             UER = len(co_eobj[f]) / len(P[len_op][op][f.rstrip('+-')])
    #             if UER < min_ei:
    #                 print('pruned')
    #                 return []

    # 搜索cesp中op各特征的演化参与实例
    for f in cesp.split(','):
        for o in co_eobj[f]:  # 对于co_obj中f的每一个实例
            # 确定o的搜索空间
            if f not in eobj or o not in eobj[f]:
                ERI = CPM_Col.searchRI(o, f, cesp, co_eobj, cgroupN, eobj)
                if ERI:
                    count_col += 1
                    for ins in ERI:
                        feature_ins = ins.split('.')[0]
                        if ins[-1] == '+' or ins[-1] == '-':
                            feature_ins = feature_ins + ins[-1]
                        if feature_ins not in eobj.keys():
                            eobj[feature_ins] = {ins}
                        else:
                            eobj[feature_ins].add(ins)
    #由于不满足反单调性，显著与否都要返回，并且已知的eobj可用于剪枝
    return eobj





# file_before = '../distances/origin_distances/test1.csv'
# file_after = '../distances/origin_distances/test2.csv'
# file_changhed = '../distances/changed_distances/test_1-2.csv'
# r = 15
# min_prev = 0.3
# min_ei = 0.4
# # points_before = get_load.read_points_from_file(file_before)
# # points_after = get_load.read_points_from_file(file_after)
# # print('points_before', points_before)
# # print('points_after', points_after)
# # P1 = CPM_Col.cpm_col(points_before, r, min_prev)
# # Biesp_Col(points_before, points_after, r, min_ei, min_prev, P1)
# instances1 = []
# instances2 = []
# changed_instances = []
# file_b = pd.read_csv(file_before)
# file_a = pd.read_csv(file_after)
# file_c = pd.read_csv(file_changhed)
# for index, row in file_b.iterrows():
#     instances1.append([row['point1'], row['point2'], row['distance']])
# for index, row in file_a.iterrows():
#     instances2.append([row['point1'], row['point2'], row['distance']])
# for index, row in file_c.iterrows():
#     changed_instances.append([row['point1'], row['point2'], row['distance']])
#
#
# points1, count1 = get_load.read_points_from_file('../testdata/data1.xlsx')
# points2, count2 = get_load.read_points_from_file('../testdata/data2.xlsx')
#
# P1 = CPM_Col.cpm_col(instances1, r, min_prev, count1)
#
#
# Biesp_Col(instances2, changed_instances, r, min_ei, min_prev, P1, count2)
