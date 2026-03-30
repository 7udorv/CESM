count_incre = 0
def m2ec_ip(instances2, changed_instances, r, min_ei, min_prev, P1, count2):

    count_candidates = 0
    P2 = copy.deepcopy(P1)
    del P2[2]
    groupN2, neighbors2 = CPM_Col.get_groupN(instances2, r)

    P2_2 = CPM_Col.gen_P2(neighbors2, count2, min_prev)

    P2[2] = P2_2

    cgroupN, neighbors = CPM_Col.get_groupN(changed_instances, r)
    esp2 = gen_esp2(neighbors)

    C2 = {}
    for esp in esp2.keys():
        p = ','.join(f.rstrip('+-') for f in esp.split(','))
        if not ('+' in esp and '-' in esp):  # 等价于 len(symbols) != 2
            if p not in C2:
                C2[p] = {}
            C2[p][esp] = esp2[esp]

    C2 = {k: v for k, v in C2.items() if k in P2_2}  # 只保留 P2_2 中存在的键的 C2
    P = C2

    k = 2
    C_P2 = CPM_Col.gen_candidates(P2_2, k + 1)

    All_Eobj = {}
    All_Eobj[k] = esp2
    next_candidates = esp2.keys()
    Sig_Esp = {}
    while True:
        cesps = gen_cesp(list(next_candidates), k, All_Eobj, P1, min_ei)
        count_candidates += len(cesps)
        next_candidates = copy.deepcopy(cesps)
        ESPs_now = {}

        for cesp in cesps:

            if '+' in cesp or '-' in cesp:  #包含变化特征
                flag, op, ep = is_OpPrevalent(cesp, P1)
                len_op = len(op.split(','))
                if flag:
                    eobj = find_eobj(cesp, k, All_Eobj, P1, op, min_ei, cgroupN)
                    if eobj == {}: 
                        next_candidates.remove(cesp)
                    if len(eobj) != 0:  #eobj不为空
                        ESPs_now[cesp] = eobj
                        for f in cesp.split(','):
                            flag_ei = True
                            if f.rstrip('+-') in op.split(','):
                                ER = len(eobj[f]) / len(P1[len_op][op][f.rstrip('+-')])
                                # op_dict[f.rstrip('+-')] = eobj[f]
                                if ER < min_ei:
                                    flag_ei = False
                                    break

                        if flag_ei:
                            if cesp not in Sig_Esp:
                                Sig_Esp[cesp] = copy.deepcopy(eobj)
                                Sig_Esp[cesp]['sequence'] = f'{op}→{ep}'
                else:  
                    next_candidates.remove(cesp)

            else:
                if len(cesp.split(',')) in P1 and cesp in P1[len(cesp.split(','))]:
                    flag = True
                else:
                    flag = False
                    next_candidates.remove(cesp)
                if flag:
                    eobj = find_eobj(cesp, k, All_Eobj, P1, cesp, min_ei, cgroupN)
                    if eobj: 
                        ESPs_now[cesp] = eobj

        k = k + 1
        All_Eobj[k] = ESPs_now
        C = {}
        if k > 2:
            # 生成候选
            C = get_update_candidates(P, C_P2, All_Eobj, cgroupN, k)
            # print('C',C)
            if not C:
                break  
            P2 = update_P2(C, C_P2, P2, k, groupN2, count2, min_prev)

        if k in P2:
            for co in list(P2[k].keys()):
                if co in C_P2:  
                    if co in C.keys():
                        for f in co.split(','):
                            obj = len(P2[k][co][f])
                            PR = obj / count2[f]
                            if PR < min_prev:  
                                del C[co]
                                del P2[k][co]
                                break
                else:
                    if co in C:
                        del C[co]
                    if co in P2[k]:
                        del P2[k][co]

            C_P2 = CPM_Col.gen_candidates(P2[k], k + 1)
            if not C_P2:
                break

            P = C  

            for cesp in list(next_candidates):
                cesp_list = cesp.split(',')
                op, ep = [], []
                for f in cesp_list:
                    if f[-1] == '+':
                        ep.append(f.rstrip('+'))
                    elif f[-1] == '-':
                        op.append(f.rstrip('-'))
                    else:
                        op.append(f)
                        ep.append(f)
                ep = ','.join(ep)
                if len(ep.split(',')) not in P2 or ep not in P2[len(ep.split(','))]:
                    if cesp in Sig_Esp.keys():
                        del Sig_Esp[cesp]
                    next_candidates.remove(cesp)
        else:  
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
                if esp1[-1].rstrip('+-') != esp2[-1].rstrip('+-'):  
                    if esp1[-1] < esp2[-1]:
                        c = esp1 + [esp2[-1]]
                    else:
                        c = esp2 + [esp1[-1]]
                    c = ','.join(c)
                    op, ep = [], []
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
    flag = False
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
    co_eobj = {}
    eobj = {}
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = All_Eobj[k][sub][f]
                co_eobj[f] = co_eobj[f] & All_Eobj[k][sub][f]

    len_op = len(op.split(','))

    if len(co_eobj) != len(cesp.split(',')): 
        return {}
    op_list = op.split(',')
    if '+' in cesp or '-' in cesp:
        for f in cesp.split(','):
            if f.rstrip('+-') in op_list:
                UER = len(co_eobj[f]) / len(P[len_op][op][f.rstrip('+-')])
                if UER < min_ei:
                    return []

    for f in cesp.split(','):
        for o in co_eobj[f]:  
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

    return eobj


def find_eobj_incre(cesp, k, All_Eobj, cgroupN):
    co_eobj = {}
    eobj = {}  
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = All_Eobj[k][sub][f]
                co_eobj[f] = co_eobj[f] & All_Eobj[k][sub][f]

    if len(co_eobj) != len(cesp.split(',')):  
        return {}


    for f in cesp.split(','):
        for o in co_eobj[f]: 
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
    return eobj


def searchRI_incre(o, f, esp, co_obj, groupN, obj):
    OSS = get_OSS_incre(o, f, esp, co_obj, groupN, obj)
    features = esp.split(',')
    features.remove(f)
    S = []
    i = 0 
    S = CPM_Col.BacktrackingSearch(S, features, OSS, i, groupN)

    if S:
        S.append(o)
    return S


def get_OSS_incre(o, f, esp, co_obj, groupN, obj): 
    features = esp.split(',')
    features.remove(f)
    OSS = {feature: [] for feature in features}
    for feature in features:
        if feature in groupN[o].keys():
            if feature in obj:
                if feature.rstrip('+-') < f:  
                    OSS[feature] = list(obj[feature] & groupN[o][feature])
                elif feature.rstrip('+-') > f:
                    if feature in groupN[o].keys():
                        if feature[-1] != '-':  
                            common_instances = co_obj[feature] & groupN[o][feature]
                        else:  
                            common_instances = groupN[o][feature]
                        obj_instances = obj[feature]
                        ordered = [x for x in common_instances if x not in obj_instances]
                        ordered.extend([x for x in common_instances if x in obj_instances])

                        OSS[feature] = ordered


            else: 
                if feature[-1] != '-':  
                    OSS[feature] = list(co_obj[feature] & groupN[o][feature])
                else:  
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

            if c1[0] < c2[0]: 
                break
            elif c1[0:k - 2] == c2[0:k - 2]:
                if c1[-1] < c2[-1]:
                    c = c1 + [c2[-1]]
                else:
                    c = c2 + [c1[-1]]
                c = ','.join(c)

                if c in C_P2:
                    esp1s_list = P[P_list[i]]
                    esp2s_list = P[P_list[j]]

                    for esp1 in esp1s_list:
                        esp1_list = esp1.split(',')
                        for esp2 in esp2s_list:
                            esp2_list = esp2.split(',')

                            if esp1_list[0:k - 2] == esp2_list[0:k - 2]:
                                if esp1_list[-1] < esp2_list[-1]:
                                    esp = esp1_list + [esp2_list[-1]]
                                else:
                                    esp = esp2_list + [esp1_list[-1]]

                                symbols = {s for x in esp for s in ['+', '-'] if s in x}
                                if len(symbols) != 2:
                                    esp = ','.join(esp)

                                    if esp in All_Eobj[k].keys(): 
                                        Eobj = copy.deepcopy(All_Eobj[k][esp])
                                    else:  
                                        Eobj = find_eobj_incre(esp, k - 1, All_Eobj, cgroupN)

                                    if Eobj:  
                                        if c not in C:
                                            C[c] = {}
                                        C[c][esp] = Eobj
    return C

def update_P2(C, C_P2, P2, k, groupN2, count2, min_prev):

    for ckey in list(C.keys()):
        if ckey in C_P2:
            if k in P2 and ckey in P2[k].keys():  
                for esp in C[ckey].keys():
                    symbols = set(ch for ch in esp if ch in '+-')
                    if symbols == {'+'}:
                        for f in esp.split(','):
                            if f.endswith('+'):
                                eobj = {item.rstrip('+') for item in C[ckey][esp][f]}
                            else:
                                eobj = C[ckey][esp][f]

                            P2[k][ckey][f.rstrip('+')] |= eobj

                    elif symbols == {'-'}:
                        items = [x.strip() for x in esp.split(',')]
                        without_minus = [x for x in items if not x.endswith('-')]
                        still_obj = {}

                        if len(without_minus) > 1:  
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
                            else:  
                                P2[k][ckey][feature] -= (C[ckey][esp][feature] - still_obj[feature])

            else:  
                obj = CPM_Col.is_prevalent(ckey, P2[k - 1], k, count2, min_prev, groupN2)

                if len(obj) == 0:  
                    del C[ckey]
                else:
                    obj = dict(sorted(obj.items()))
                    if obj:
                        if k not in P2:
                            P2[k] = {}
                        P2[k][ckey] = obj

    return P2
