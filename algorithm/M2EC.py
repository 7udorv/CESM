count_col = 0

def m2ec(instances2, changed_instances, r, min_ei, min_prev, P1, count2):

    count_candidates = 0
    P2 = CPM_Col.cpm_col(instances2, r, min_prev, count2)
    cgroupN, neighbors = CPM_Col.get_groupN(changed_instances, r)
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
            if '+' in cesp or '-' in cesp: 
                if flag:

                    eobj = find_eobj(cesp, k, All_Eobj, P1, op, min_ei, cgroupN)
                    # print('eobj', eobj)
                    if eobj == {}: 
                        next_candidates.remove(cesp)
                    if len(eobj) != 0:  
                        ESPs_now[cesp] = eobj
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
                        if flag_ei:  
                            if cesp not in Sig_Esp:
                                Sig_Esp[cesp] = eobj
                                Sig_Esp[cesp]['sequence'] = f'{op}→{ep}'
                                # Sig_Esp[cesp]['op_dict'] = op_dict
                                # Sig_Esp[cesp]['ep_dict'] = ep_dict

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
                if esp1[-1].rstrip('+-') != esp2[-1].rstrip('+-'): 
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
    global count_col
    co_eobj = {}
    eobj = {}  
    subsets = [','.join(sub) for sub in combinations(cesp.split(','), k)]
    for sub in subsets:
        if sub in All_Eobj[k].keys():
            features = sub.split(',')
            for f in features:
                if f not in co_eobj.keys():
                    co_eobj[f] = set()
                co_eobj[f] = co_eobj[f] | All_Eobj[k][sub][f]
    len_op = len(op.split(','))

    if len(co_eobj) != len(cesp.split(',')): 
        return {}
    for f in cesp.split(','):
        for o in co_eobj[f]:  
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

    return eobj


