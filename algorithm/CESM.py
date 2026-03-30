def CESM(count_dict, original_distances_dict, changed_distances_dict, r, min_prev, min_ei, min_overlap, original_files):
    i = 0
    ESPs = {}
    count_list = list(count_dict.keys())
    o_distances_list = list(original_distances_dict.keys())
    c_distances_list = list(changed_distances_dict.keys())
    while True:
        print(i)
        if i == 0:
            instances1 = original_distances_dict[o_distances_list[i]]
            count1 = count_dict[count_list[i]]
            P1 = CPM_Col.cpm_col(instances1, r, min_prev, count1)
            # print('P1',P1)

        instances2 = original_distances_dict[o_distances_list[i+1]]
        count2 = count_dict[count_list[i+1]]

        changed_instances = changed_distances_dict[c_distances_list[i]]

        Sig_ESP, P2 = M2EC.m2ec(instances2, changed_instances, r, min_ei, min_prev, P1, count2)
        ESPs[i] = Sig_ESP
        P1 = P2
        count1 = count2

        if i+1 >= len(original_files)-1:
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
                        op_dict = {}
                        ep_dict = {}
                        for f1 in c1_ep.split(','):
                            for f2 in esp1.split(','):
                                if f1 == f2.rstrip("+-"):
                                    ins_ep = ESPs[n][esp1][f2]
                                    if f1 == f2:  
                                        ep_dict[f1] = ins_ep
                                    else:
                                        ins_ep = set(ins.rstrip("+-") for ins in ins_ep)
                                        ep_dict[f1] = ins_ep
                            for f3 in esp2.split(','):
                                if f1 == f3.rstrip("+-"):
                                    ins_op = ESPs[n + 1][esp2][f3]
                                    if f1 == f3:  
                                        op_dict[f1] = ins_op
                                    else:
                                        ins_op = set(ins.rstrip("+-") for ins in ins_op)
                                        op_dict[f1] = ins_op

                            IOI = len(ep_dict[f1] & op_dict[f1]) / len(ep_dict[f1])
                            if IOI < min_ioi:
                                flag_ioi = False
                                break
                        if flag_ioi:
                            parts1 = c1.split("→")
                            parts2 = c2.split("→")
                            s = "→".join(parts1 + [parts2[-1]])
                            LS.add(s)
                        else: 
                            if n not in P:
                                P[n] = []
                            P[n].append(c1)
            if not flag_connectable: 
                if n not in P:
                    P[n] = []
                P[n].append(c1)

        C1 = C2 + list(LS)
        C1 = sorted(C1)
        n = n + 1
        if n + 1 > len(original_files) - 2:
            P[n] = C1
            break
    return P


