from copy import deepcopy
from itertools import combinations
import copy
from ESPM_Col import CPM_Col_incre
import get_load
from ESPM_Col import CPM_Col
import re

count_incre = 0
def Biesp_Col_incre(instances2, changed_instances, r, min_ei, min_prev, P1, count2):

    count_candidates = 0
    #首先获取P2_2
    P2 = copy.deepcopy(P1)
    del P2[2]
    # 计算P2的groupN2，并生成neighbors
    groupN2, neighbors2 = CPM_Col.get_groupN(instances2, r)

    # 生成P2的二阶
    P2_2 = CPM_Col.gen_P2(neighbors2, count2, min_prev)

    P2[2] = P2_2

    # 获取所有变化点及其邻居的邻居
    cgroupN, neighbors = CPM_Col.get_groupN(changed_instances, r)
    # 根据neighbors生成所有二阶变化同位模式,其中包括不包含变化特征的模式
    esp2 = gen_esp2(neighbors)

    #获取变化的模式
    # 首先利用二阶变化模式生成二阶候选
    C2 = {}
    for esp in esp2.keys():
        # 只要特征相同就加入，无论是否变化
        p = ','.join(f.rstrip('+-') for f in esp.split(','))
        # symbols = set(re.findall(r'[+-]', esp))
        # if len(symbols) != 2:  # 过滤包含不同变化特征的模式
        #     if p not in C2:
        #         C2[p] = {}
        #     C2[p][esp] = esp2[esp]
        if not ('+' in esp and '-' in esp):  # 等价于 len(symbols) != 2
            if p not in C2:
                C2[p] = {}
            C2[p][esp] = esp2[esp]

    # 将不频繁的P2二阶从C2中删去
    C2 = {k: v for k, v in C2.items() if k in P2_2}  # 只保留 P2_2 中存在的键的 C2
    P = C2

    k = 2
    # 生成P2的三阶候选
    C_P2 = CPM_Col.gen_candidates(P2_2, k + 1)

    All_Eobj = {}
    All_Eobj[k] = esp2
    next_candidates = esp2.keys()
    Sig_Esp = {}
    while True:
        # 生成候选演化模式以及其演化参与实例
        # cesps = gen_cesp(list(next_candidates), k)
        cesps = gen_cesp(list(next_candidates), k, All_Eobj, P1, min_ei)
        count_candidates += len(cesps)
        next_candidates = copy.deepcopy(cesps)
        ESPs_now = {}

        for cesp in cesps:
            #首先检查cesp中是否包含变化特征，不包含则直接加入候选
            if '+' in cesp or '-' in cesp:  #包含变化特征
                #检查原模式是否频繁，对于在P中的，直接认为op频繁，对于不在P中的，需要进一步检验频繁性
                flag, op, ep = is_OpPrevalent(cesp, P1)
                len_op = len(op.split(','))
                #记录已找到的变化op
                if flag:
                    #获取cesp的所有演化参与实例
                    eobj = find_eobj(cesp, k, All_Eobj, P1, op, min_ei, cgroupN)
                    if eobj == {}:  #表示没有参与实例，删除候选
                        next_candidates.remove(cesp)
                    if len(eobj) != 0:  #eobj不为空
                        ESPs_now[cesp] = eobj
                        #检查cesp是否满足显著阈值
                        for f in cesp.split(','):
                            flag_ei = True
                            if f.rstrip('+-') in op.split(','):
                                ER = len(eobj[f]) / len(P1[len_op][op][f.rstrip('+-')])
                                # op_dict[f.rstrip('+-')] = eobj[f]
                                if ER < min_ei:
                                    flag_ei = False
                                    break
                            # if f.rstrip('+-') in ep.split(','):
                            #     ep_dict[f.rstrip('+-')] = eobj[f]
                        if flag_ei:  #保存所有op频繁，并且满足min_ei的cesp
                            #那么cesp是显著的
                            if cesp not in Sig_Esp:
                                Sig_Esp[cesp] = copy.deepcopy(eobj)
                                Sig_Esp[cesp]['sequence'] = f'{op}→{ep}'
                                # Sig_Esp[cesp]['op_dict'] = op_dict
                                # Sig_Esp[cesp]['ep_dict'] = ep_dict
                else:  #如果op不频繁，从候选中删除cesp
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

        k = k + 1

        All_Eobj[k] = ESPs_now
        # 更新k时的P
        C = {}
        if k > 2:
            # 生成候选
            C = get_update_candidates(P, C_P2, All_Eobj, cgroupN, k)
            # print('C',C)
            if not C:
                break  #说明已经没有频繁ep
            # 更新ep的参与实例，对于消失实例，需要检查是否有共享，新增实例不管
            P2 = update_P2(C, C_P2, P2, k, groupN2, count2, min_prev)

        # 更新完毕，获取在候选中的，并计算
        if k in P2:
            for co in list(P2[k].keys()):
                if co in C_P2:  #只计算在候选的
                    # 判断有没有变化
                    if co in C.keys():  #变化的才计算
                        for f in co.split(','):
                            obj = len(P2[k][co][f])
                            PR = obj / count2[f]
                            if PR < min_prev:  #如果不频繁
                                #从P2中删除co
                                del C[co]
                                del P2[k][co]
                                break
                else:
                    # 不在候选则删除
                    if co in C:
                        del C[co]
                    if co in P2[k]:
                        del P2[k][co]

            #生成k+1阶候选
            C_P2 = CPM_Col.gen_candidates(P2[k], k + 1)
            if not C_P2:
                break

            P = C  # P用于生成变化候选
            # print(P2)

            #根据Sig_Esp中所有模式获取P2
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
        else:  #k不在P2中
            #删除所有k对应的sig，并break
            for cesp in list(next_candidates):
                if cesp in Sig_Esp.keys():
                    del Sig_Esp[cesp]
                    next_candidates.remove(cesp)
            break
        if len(next_candidates) == 0:
            break
    print(count_incre)
    print('count_candidates', count_candidates)
    return Sig_Esp, P2


def gen_esp2(neighbors):
    candidates = {}
    for neighbor in neighbors:
        pt1 = neighbor[0]
        f1 = pt1.split('.')[0]
        pt2 = neighbor[1]
        f2 = pt2.split('.')[0]
        if pt1[-1] == '+' or pt1[-1] == '-':
            f1 = f1 + pt1[-1]
        if pt2[-1] == '+' or pt2[-1] == '-':
            f2 = f2 + pt2[-1]
        features = f1 + ',' + f2
        if features not in candidates.keys():
            candidates[features] = {}
            candidates[features][f1] = set()
            candidates[features][f2] = set()
        candidates[features][f1].add(pt1)
        candidates[features][f2].add(pt2)
    return candidates


# def gen_cesp(esp_list, k):
#     esp_list = sorted(esp_list, key=lambda x: x[0])
#     cesp = set()
#     for i in range(0, len(esp_list) - 1):
#         esp1 = esp_list[i].split(',')
#         for j in range(i + 1, len(esp_list)):
#             esp2 = esp_list[j].split(',')
#             if esp1[0].rstrip('+-') < esp2[0].rstrip('+-'):
#                 break
#             elif esp1[0:k - 1] == esp2[0:k - 1]:
#                 if esp1[-1].rstrip('+-') != esp2[-1].rstrip('+-'):  #要求最后一个元素的特征不能相同
#                     if esp1[-1] < esp2[-1]:
#                         c = esp1 + [esp2[-1]]
#                     else:
#                         c = esp2 + [esp1[-1]]
#                     c = ','.join(c)
#                     cesp.add(c)
#     return cesp

def gen_cesp(esp_list, k, All_Eobj, P1, min_ei):
    esp_list = sorted(esp_list, key=lambda x: x[0])
    cesp = set()
    for i in range(0, len(esp_list) - 1):
        for j in range(i + 1, len(esp_list)):
            esp1 = esp_list[i].split(',')
            esp2 = esp_list[j].split(',')
            if esp1[0].rstrip('+-') < esp2[0].rstrip('+-'):
                break
            elif esp1[0:k - 1] == esp2[0:k - 1]:
                if esp1[-1].rstrip('+-') != esp2[-1].rstrip('+-'):  #要求最后一个元素的特征不能相同
                    if esp1[-1] < esp2[-1]:
                        c = esp1 + [esp2[-1]]
                    else:
                        c = esp2 + [esp1[-1]]
                    c = ','.join(c)
                    #判断子集的频繁性，key是子集
                    op, ep = [], []
                    # 首先寻找op和ep
                    for f in c.split(','):
                        if f[-1] == '+':
                            ep.append(f.rstrip('+'))
                        elif f[-1] == '-':
                            op.append(f.rstrip('-'))
                        else:
                            op.append(f)
                            ep.append(f)
                    op = ','.join(op)
                    flag_s_1 = is_SPR(esp_list[i], All_Eobj, P1, op, k, min_ei)
                    if flag_s_1:
                        flag_s_2 = is_SPR(esp_list[j], All_Eobj, P1, op, k, min_ei)
                        if flag_s_1 and flag_s_2:
                            cesp.add(c)
    return cesp


def is_SPR(key, All_Eobj, P1, op, k, min_ei):
    flag_suer = True
    if len(op.split(',')) > 1 and k in All_Eobj and key in All_Eobj[k]:
        for f in key.split(','):
            if len(op.split(',')) in P1:
                if op in P1[len(op.split(','))] and f.rstrip('+-') in P1[len(op.split(','))][op]:
                    super_UER = len(All_Eobj[k][key][f]) / len(P1[len(op.split(','))][op][f.rstrip('+-')])
                    if super_UER < min_ei:
                        # print("super_pruned")
                        flag_suer = False
                        return flag_suer

    return flag_suer


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
    op = ','.join(op)
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
    global count_incre
    #获取cesp的候选演化参与实例
    # 获取候选参与实例
    co_eobj = {}
    eobj = {}  #用于存放所有演化模式的演化参与实例
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = All_Eobj[k][sub][f]
                co_eobj[f] = co_eobj[f] & All_Eobj[k][sub][f]

    #获取op的obj
    len_op = len(op.split(','))

    if len(co_eobj) != len(cesp.split(',')):  #说明缺失的特征没有演化参与实例
        return {}
    op_list = op.split(',')
    # 计算演化参与率上界(只有演化同位模式需要)
    if '+' in cesp or '-' in cesp:
        for f in cesp.split(','):
            if f.rstrip('+-') in op_list:
                UER = len(co_eobj[f]) / len(P[len_op][op][f.rstrip('+-')])
                if UER < min_ei:
                    # print('pruned')
                    return []

    # 搜索cesp中op各特征的演化参与实例
    for f in cesp.split(','):
        for o in co_eobj[f]:  # 对于co_obj中f的每一个实例
            # 确定o的搜索空间
            if f not in eobj or o not in eobj[f]:
                ERI = CPM_Col.searchRI(o, f, cesp, co_eobj, cgroupN, eobj)
                if ERI:
                    count_incre += 1
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


def find_eobj_incre(cesp, k, All_Eobj, cgroupN):
    #获取cesp的候选演化参与实例
    # 获取候选参与实例
    co_eobj = {}
    eobj = {}  #用于存放所有演化模式的演化参与实例
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = All_Eobj[k][sub][f]
                co_eobj[f] = co_eobj[f] & All_Eobj[k][sub][f]

    if len(co_eobj) != len(cesp.split(',')):  #说明缺失的特征没有演化参与实例
        return {}

    # 搜索cesp中op各特征的演化参与实例
    for f in cesp.split(','):
        for o in co_eobj[f]:  # 对于co_obj中f的每一个实例
            # 确定o的搜索空间
            if f not in eobj or o not in eobj[f]:
                ERI = CPM_Col.searchRI(o, f, cesp, co_eobj, cgroupN, eobj)
                if ERI:
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


def searchRI_incre(o, f, esp, co_obj, groupN, obj):
    # 获取搜索空间
    OSS = get_OSS_incre(o, f, esp, co_obj, groupN, obj)
    features = esp.split(',')
    features.remove(f)
    S = []
    i = 0  # 用于记录特征位置
    S = CPM_Col.BacktrackingSearch(S, features, OSS, i, groupN)

    # 如果S不为空
    if S:
        S.append(o)
    return S


def get_OSS_incre(o, f, esp, co_obj, groupN, obj):  #只用搜索未变化特征
    features = esp.split(',')
    features.remove(f)
    # 确认每一个特征的搜索空间
    OSS = {feature: [] for feature in features}
    for feature in features:
        if feature in groupN[o].keys():
            # 如果obj已知
            if feature in obj:
                if feature.rstrip('+-') < f:  # 对于一个ABC,如果此时在查找B，A的的所有参与实例一定已经搜索完全，此时obj已知
                    OSS[feature] = list(obj[feature] & groupN[o][feature])
                elif feature.rstrip('+-') > f:  # feature 中搜索到了部分obj，也算作obj未知
                    # 把未搜索的部分放在最后
                    if feature in groupN[o].keys():
                        if feature[-1] != '-':  #对于未变化特征
                            # 求交集
                            common_instances = co_obj[feature] & groupN[o][feature]
                        else:  #变化特征没有co_obj
                            common_instances = groupN[o][feature]
                        # 将 obj[feature] 中的实例移动到末尾
                        obj_instances = obj[feature]
                        # 保持顺序的同时将 obj_instances 中的先剔除再添加到末尾
                        ordered = [x for x in common_instances if x not in obj_instances]
                        ordered.extend([x for x in common_instances if x in obj_instances])

                        # 更新 OSS[feature]
                        OSS[feature] = ordered


            else:  # obj完全未知
                if feature[-1] != '-':  #未变化特征
                    OSS[feature] = list(co_obj[feature] & groupN[o][feature])
                else:  # 变化特征
                    OSS[feature] = list(groupN[o][feature])
        else:
            OSS[feature] = []
    return OSS

def get_update_candidates(P, C_P2, All_Eobj, cgroupN, k):

    C = {}
    P_list = list(P.keys())

    for i in range(0, len(P_list) - 1):
        c1 = P_list[i].split(',')
        for j in range(i + 1, len(P_list)):
            c2 = P_list[j].split(',')

            if c1[0] < c2[0]:  # A<B 直接跳过
                break
            # 如果k-1相同就连接
            elif c1[0:k - 2] == c2[0:k - 2]:
                # 比较最后一个元素大小并连接
                if c1[-1] < c2[-1]:
                    c = c1 + [c2[-1]]
                else:
                    c = c2 + [c1[-1]]
                c = ','.join(c)

                # 检查候选是否在C_P2
                if c in C_P2:
                    esp1s_list = P[P_list[i]]
                    esp2s_list = P[P_list[j]]

                    for esp1 in esp1s_list:
                        esp1_list = esp1.split(',')
                        for esp2 in esp2s_list:
                            esp2_list = esp2.split(',')

                            # 如果k-1相同就连接
                            if esp1_list[0:k - 2] == esp2_list[0:k - 2]:
                                # 比较最后一个元素大小并连接
                                if esp1_list[-1] < esp2_list[-1]:
                                    esp = esp1_list + [esp2_list[-1]]
                                else:
                                    esp = esp2_list + [esp1_list[-1]]

                                # 检查是否包含同时 + 和 -
                                symbols = {s for x in esp for s in ['+', '-'] if s in x}
                                if len(symbols) != 2:
                                    esp = ','.join(esp)

                                    if esp in All_Eobj[k].keys():  # 已存在
                                        Eobj = copy.deepcopy(All_Eobj[k][esp])
                                    else:  # 需要增量计算
                                        Eobj = find_eobj_incre(esp, k - 1, All_Eobj, cgroupN)

                                    if Eobj:  # 非空才保存
                                        if c not in C:
                                            C[c] = {}
                                        C[c][esp] = Eobj
    return C

def update_P2(C, C_P2, P2, k, groupN2, count2, min_prev):

    for ckey in list(C.keys()):  # C2 中的都是变化的，都要更新
        if ckey in C_P2:  # 如果在 k 阶候选中
            if k in P2 and ckey in P2[k].keys():  # 如果 ckey 在 P2[k] 中
                for esp in C[ckey].keys():  # 对于 ckey 对应的每一个 esp
                    # 判断 esp 的符号
                    symbols = set(ch for ch in esp if ch in '+-')

                    # ----------- 新增 ----------
                    if symbols == {'+'}:
                        for f in esp.split(','):
                            if f.endswith('+'):
                                eobj = {item.rstrip('+') for item in C[ckey][esp][f]}
                            else:
                                eobj = C[ckey][esp][f]

                            P2[k][ckey][f.rstrip('+')] |= eobj

                    # ----------- 删除 ----------
                    elif symbols == {'-'}:
                        items = [x.strip() for x in esp.split(',')]
                        without_minus = [x for x in items if not x.endswith('-')]
                        still_obj = {}

                        if len(without_minus) > 1:  # 还剩至少两个未删除特征才能形成边
                            for f in without_minus:
                                for o in C[ckey][esp][f]:
                                    if f not in still_obj or o in still_obj[f]:
                                        RI = searchRI_incre(o, f, esp, C[ckey][esp], groupN2, still_obj)
                                        if RI:
                                            for ins in RI:
                                                feature_ins = ins.split('.')[0]
                                                if ins[-1] in ['+', '-']:
                                                    feature_ins = feature_ins + ins[-1]
                                                if feature_ins not in still_obj:
                                                    still_obj[feature_ins] = {ins}
                                                else:
                                                    still_obj[feature_ins].add(ins)

                        for feature in esp.split(','):
                            if feature.endswith('-') or still_obj == {}:
                                cleaned_objs = {obj.rstrip('-') for obj in C[ckey][esp][feature]}
                                P2[k][ckey][feature.rstrip('-')] -= cleaned_objs
                            else:  # 未变化
                                P2[k][ckey][feature] -= (C[ckey][esp][feature] - still_obj[feature])

            else:  # 直接从 k-1 层中找
                obj = CPM_Col.is_prevalent(ckey, P2[k - 1], k, count2, min_prev, groupN2)

                if len(obj) == 0:  # 不频繁则删除该键
                    del C[ckey]
                else:
                    obj = dict(sorted(obj.items()))
                    if obj:
                        if k not in P2:
                            P2[k] = {}
                        P2[k][ckey] = obj

    return P2

# file_before = '../testdata/data1.xlsx'
# file_after = '../testdata/data2.xlsx'
# r = 15
# min_prev = 0.3
# min_ei = 0.4
# points_before = get_load.read_points_from_file(file_before)
# points_after = get_load.read_points_from_file(file_after)
# print('points_before', points_before)
# print('points_after', points_after)
# P1 = CPM_Col.cpm_col(points_before, r, min_prev)
# Biesp_Col_incre(points_before, points_after, r, min_ei, min_prev, P1)
