import desc as dsc
import splitClass as sc
import statTests as stat
from math import fabs
from math import sqrt
from scipy.stats import norm
import networkx as nx


root = '0'


def subgroup_identification(sample, a_logrank, cov_at_level):
    subgroups = nx.DiGraph()
    cov_used = {key: val for key in sample.keys() for val in [False]}
    subgroups.add_node(root)
    sub_ident(sample, subgroups, cov_used, a_logrank, cov_at_level, parent=root, recurs_level=1)
    return subgroups


def get_group(groups, node_name):
    return groups.node[node_name][dsc.node_content]


def sub_ident(sample, subgroups, cov_used, a_logrank, cov_at_level, parent, recurs_level):
    def sort_by_p_val(itm):
        return itm[-1]

    def merge(grp, grp0):
        grp0.extend(grp)
        grp0.sort(key=sort_by_p_val)
        k = len(grp0)
        while k > cov_at_level:
            grp0.pop()
            k -= 1
        return

    def create_node_name(rec_lvl, num_nd, mode):
        if not (mode == 'left' or mode == 'right'):
            exit('creat_node_name() error')
        return '-'.join([str(rec_lvl), str(num_nd), mode])

    def better_group(grp1, grp2):
        return 1 if grp1.logrank.p_value < grp2.logrank.p_value else -1

    def select_subsample(samp, group_inds):
        subsamp = {key: val for key in samp.keys() for val in [[]]}
        for cov in subsamp.keys():
            subsamp[cov] = [samp[cov][i] for i in xrange(len(samp[cov])) if i in group_inds]
        return subsamp

    group0 = []
    for cov_name, val in cov_used.items():
        if val or cov_name == dsc.covariates2002[0]\
               or cov_name == dsc.covariates2002[1]\
               or cov_name == dsc.covariates2002[2]:
            continue
        group = select_best_k_splits(sample, cov_name, a_logrank, cov_at_level)
        merge(group, group0)
    gamma = select_cont_param()
    better = []
    groups_next = []
    cov_next = []
    for i, split in enumerate(group0):
        name1, name2 = create_node_name(recurs_level, i, 'left'), create_node_name(recurs_level, i, 'right')
        subgroups.add_node(name1, {dsc.node_content: split[0]})
        subgroups.add_edge(parent, name1, {dsc.edge_content: split[1]})
        subgroups.add_node(name2, {dsc.node_content: split[2]})
        subgroups.add_edge(parent, name2, {dsc.edge_content: split[3]})
        better.append(
            name1 if better_group(get_group(subgroups, name1), get_group(subgroups, name2)) > 0 else name2
        )
        if parent == root:
            groups_next.append(better[i])
            cov_next.append(split[4])
        else:
            if continuation_criterion(get_group(subgroups, better[i]), get_group(subgroups, parent), gamma):
                groups_next.append(better[i])
                cov_next.append(split[4])
    for i in xrange(len(groups_next)):
        subsample = select_subsample(sample, get_group(subgroups, groups_next[i]).group)
        new_cov_used = cov_used.copy()
        new_cov_used[cov_next[i]] = True
        sub_ident(subsample, subgroups, new_cov_used, a_logrank, cov_at_level, [i], recurs_level + 1)
    return subgroups


def continuation_criterion(current, parent, gamma):
    return current.logrank.p_value <= gamma * parent.logrank.p_value


def select_cont_param():
    return 1


def select_best_k_splits(sample, cov_name, a_logrank, cov_at_level):
    # returning values:
    # group = [ tuple(group1, desc1, group2, desc2, cov_name, p_value_split), ... ]

    def get_cov_values(samp, cov):
        return tuple(set(samp[cov]))

#    def get_num_splits(cov_name, values):
#        if cov_name == 'Sex' or cov_name == 'Mediastinum' or cov_name == 'CNS':
#            return 1
#        elif cov_name == 'Age' or cov_name == 'Leuc' or cov_name == 'Leber' or cov_name == 'Milz':
#            return len(values) - 1
#        elif cov_name == 'Immun':
#            return 2 * pow(2, len(values) / 2)
#        else:
#            exit('get_num_splits() error')

    def make_split(samp, cov, vals):
        if cov == 'Sex':
            l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                   [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
            yield l, r, 'Male', 'Female'
        elif cov == 'Immun':
            immun_b, immun_t = vals[:len(vals) / 2], vals[len(vals) / 2:]
            for k in xrange(1, 2 ** len(immun_b)):
                l, r = [], []
                s = bin(k)[2:]
                if len(s) < len(immun_b):
                    s = "".join(['0' for j in xrange(len(immun_b) - len(s))]) + s
                level = tuple([immun_b[j] for j in xrange(len(s)) if s[j] == '1'])
                not_level = tuple([immun_b[j] for j in xrange(len(s)) if s[j] == '0'])
                for j in xrange(len(samp[cov])):
                    if samp[cov][j] in level:
                        r.append(j)
                    else:
                        l.append(j)
                yield l, r, 'in %s' % (str(level)), 'in %s' % (str(not_level))
            for k in xrange(1, 2 ** len(immun_t)):
                l, r = [], []
                s = bin(k)[2:]
                if len(s) < len(immun_t):
                    s = "".join(['0' for j in xrange(len(immun_t) - len(s))]) + s
                level = tuple([immun_t[j] for j in xrange(len(s)) if s[j] == '1'])
                not_level = tuple([immun_t[j] for j in xrange(len(s)) if s[j] == '0'])
                for j in xrange(len(samp[cov])):
                    if samp[cov][j] in level:
                        r.append(j)
                    else:
                        l.append(j)
                yield l, r, 'Immun = %s' % (str(level)), 'Immun = %s' % (str(not_level))
        elif cov == 'CNS':
            l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                          [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
            yield l, r, 'CNS = 1', 'CNS = 2'
        elif cov == 'Mediastinum':
            l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                          [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
            yield l, r, 'Mediastinum = 1', 'Mediastinum = 2'
        elif cov == 'Age' or cov == 'Leuc' or cov == 'Leber' or cov == 'Milz':
            for level in vals:
                l, r = [], []
                for j in xrange(len(samp[cov])):
                    if samp[cov][j] < level:
                        l.append(j)
                    else:
                        r.append(j)
                yield l, r, '%s >= %f' % (cov, level), '%s < %f' % (cov, level)
        else:
            exit('make_split() error')

    def convert_data(outcm):
        # 0 - patient is alive, 1 - data lost
        censor = [j for j in xrange(len(outcm)) if outcm[j] == 0 or outcm[j] == 1]
        convert_outcm = [True for x in outcm]
        for j in censor:
            convert_outcm[j] = False
        return convert_outcm

    def sep_treats(outcm, tm, trt, trt_type):
        return [outcm[j] for j in xrange(len(outcm)) if trt[j] == trt_type[0]],\
               [tm[j] for j in xrange(len(outcm)) if trt[j] == trt_type[0]],\
               [outcm[j] for j in xrange(len(outcm)) if trt[j] == trt_type[1]],\
               [tm[j] for j in xrange(len(outcm)) if trt[j] == trt_type[1]]

    def add_split(splts, other_splt):
        pos = -1
        j = 0
        while j < len(splts):
            if splts[j][-1] > other_splt[-1]:
                pos = j
                break
            j += 1
        if pos == -1:
            splts.append(other_splt)
        else:
            splts.insert(pos, other_splt)
        return splts

#    def adjust_pval(p_split, num_splits):
#        return 1 - pow(1 - p_split, num_splits)

    treatment = sample[dsc.covariates2002[0]]  # Rand - protocol
    outcome = sample[dsc.covariates2002[1]]  # Tod - alive, lost, dead
    time = sample[dsc.covariates2002[2]]  # Time - lifetime
    treatment_type = sorted(list(set(treatment)))
    best_k_splits = []
    cov_values = get_cov_values(sample, cov_name)
    if len(cov_values) == 1:
        return best_k_splits
    for left, right, desc_left, desc_right in make_split(sample, cov_name, cov_values):
        if len(left) < dsc.min_sub_size or len(right) < dsc.min_sub_size:
            continue
        group1, group2 = sc.Subgroup(), sc.Subgroup()
        treat1, treat2 = [treatment[i] for i in left], [treatment[i] for i in right]
        out1, out2 = [outcome[i] for i in left], [outcome[i] for i in right]
        time1, time2 = [time[i] for i in left], [time[i] for i in right]
        convert_out1, convert_out2 = convert_data(out1), convert_data(out2)
        out11, t11, out12, t12 = sep_treats(convert_out1, time1, treat1, treatment_type)
        out21, t21, out22, t22 = sep_treats(convert_out2, time2, treat2, treatment_type)
        kmf11 = stat.kaplan_meier(out11, t11, treatment_type[0])
        kmf12 = stat.kaplan_meier(out12, t12, treatment_type[1])
        kmf21 = stat.kaplan_meier(out21, t21, treatment_type[0])
        kmf22 = stat.kaplan_meier(out22, t22, treatment_type[1])
        surv11 = stat.get_kmf_survival(kmf11)
        surv12 = stat.get_kmf_survival(kmf12)
        surv21 = stat.get_kmf_survival(kmf21)
        surv22 = stat.get_kmf_survival(kmf22)
        res1 = stat.logrank(out11, t11, out12, t12, alpha=a_logrank)
        res2 = stat.logrank(out21, t21, out22, t22, alpha=a_logrank)
        pow1 = stat.logrank_power(min(len(out11), len(out12)), surv11, surv12, alpha=a_logrank)
        pow2 = stat.logrank_power(min(len(out21), len(out22)), surv21, surv22, alpha=a_logrank)
        p_split = split_criterion(res1.test_statistic, res2.test_statistic)
        group1.set_subgroup(left, kmf11, treatment_type[0], kmf12, treatment_type[1], logrank=res1, pwr=pow1)
        group2.set_subgroup(right, kmf21, treatment_type[0], kmf22, treatment_type[1], logrank=res2, pwr=pow2)
        add_split(best_k_splits, (group1, desc_left, group2, desc_right, cov_name, p_split))
        if len(best_k_splits) > cov_at_level:
            best_k_splits.pop()
    return best_k_splits


def split_criterion(zstat1, zstat2, mode=1):
    def crit1(z1, z2):
        return 2 * (1 - norm.cdf(fabs(z1 - z2) / sqrt(2)))

    def crit2(z1, z2):
        return 2 * min(1 - norm.cdf(z1), 1 - norm.cdf(z2))

    def crit3(z1, z2):
        return max(crit1(z1, z2), crit2(z1, z2))

    if mode == 1:
        f = crit1
    elif mode == 2:
        f = crit2
    else:
        f = crit3
    return f(zstat1, zstat2)


"""
def sub_ident(s, inds, covs, treat, outs, time, tree_h, recurs_level, parent_subgr=None):

    def max_num_splits():
        return dsc.numBestCov

    def set_subgroup(sp_subgroup, sp_z_val, sp_p_logrank, sp_pwr, sp_bet, sp_surv, sp_n_surv, sp_n_total,
                     sp_p_crit, sp_lev, ind, covs):
        gr = spCl.Subgroup()

        gr.group = sp_subgroup
        gr.cov_ind = ind
        gr.cov_name = covs[ind][0]
        gr.cov_level = sp_lev
        gr.z = sp_z_val
        gr.p_val_logrank = sp_p_logrank
        gr.pwr = sp_pwr
        gr.p_val_crit = sp_p_crit
        gr.sz = sp_n_total
        gr.surv_n = sp_n_surv
        gr.surv = sp_surv
        gr.better = sp_bet

        return gr

    def make_subgroups(sp_subgroup, sp_z_val, sp_p_logrank, sp_pwr, sp_bet, sp_surv, sp_n_surv, sp_n_total,
                       sp_p_crit, sp_lev, ind, covs):
        groups = []
        for i in xrange(len(sp_subgroup)):
            groups.append(set_subgroup(sp_subgroup[i], sp_z_val[i], sp_p_logrank[i], sp_pwr[i], sp_bet[i],
                                       sp_surv[i], sp_n_surv[i], sp_n_total[i], sp_p_crit[i], sp_lev[i], ind, covs))
        return groups

    def insert_into_list(dest, src, max_len=max_num_splits()):

        def merge(dest, src):
            dest.extend(src)
            dest.sort(key=spCl.sort_by_p_val)
            return

        src.sort(key=spCl.sort_by_p_val)
        merge(dest, src[:max_len + 1])
        return dest[:max_len + 1]

    def set_descr(subgr):
        if subgr.desc is None:
            subgr.desc = 'sz1: %4d, sz2: %4d, surv1: %4d, surv2: %4d, cov: %12s, lev: %20s, bet: %d, Z: %8.4f, p-val: %.4f, pow: ' % \
                (subgr.sz[0], subgr.sz[1], subgr.surv_n[0], subgr.surv_n[1], subgr.cov_name, str(subgr.cov_level), subgr.better, subgr.z, subgr.p_val_logrank)
        else:
            subgr.desc += '\nsz1: %4d, sz2: %4d, surv1: %4d, surv2: %4d, cov: %12s, lev: %20s, bet: %d, Z: %8.4f, p-val: %.4f, pow: ' % \
                (subgr.sz[0], subgr.sz[1], subgr.surv_n[0], subgr.surv_n[1], subgr.cov_name, str(subgr.cov_level), subgr.better, subgr.z, subgr.p_val_logrank)
        subgr.desc += 'None' if subgr.pwr is None else '%.4f' % (subgr.pwr)
        return subgr

    def make_descs(groups):
        for i in xrange(len(groups)):
            set_descr(groups[i])
        return groups

    #res_sub, res_z, res_p_logrank, res_pwr, res_bet, res_surv = [], [], [], [], [], []
    #res_s_crit = []
    prom_subgroups = []
    res_groups = []

    print 'recurs_level: %d' % (recurs_level)

    if recurs_level < tree_h:
        for cov in inds:

            print 'cov: %d ' % (cov)

            sp_subgroup, sp_z_val, sp_p_logrank, sp_pwr, sp_bet, sp_surv, sp_n_surv, sp_n_total, sp_p_crit, sp_lev = \
                select_subgroup(s, cov, covs, treat, outs, time)

            subgroups = make_subgroups(sp_subgroup, sp_z_val, sp_p_logrank, sp_pwr, sp_bet,
                                       sp_surv, sp_n_surv, sp_n_total, sp_p_crit, sp_lev, cov, covs)

            #insert_into_list(res_sub, res_z, res_p_logrank, res_pwr, res_bet, res_surv,
            #                 sp_subgroup, sp_z_val, sp_p_logrank, sp_pwr, sp_bet, sp_surv, sp_p_crit)
            prom_subgroups = insert_into_list(prom_subgroups, subgroups)

        # list of pairs: [0] subgroup, [1] either left child (1) or right one (2)
        #best = select_best(groups)

        print 'num_prom_subgroups: %d' %(len(prom_subgroups))

        prom_subgroups = make_descs(prom_subgroups)

        subgroups_to_split = continuation_criterion(prom_subgroups, parent_subgr)
        subgroups_other = set(range(len(prom_subgroups))) - set(subgroups_to_split)

        print 'groups to split: %d, groups result: %d' % (len(subgroups_to_split), len(subgroups_other))

        res_groups.extend([prom_subgroups[i] for i in subgroups_other])

        for i in subgroups_to_split:
            new_gr = prom_subgroups[i].group
            new_s = [s[j] for j in new_gr]
            new_treat = [treat[j] for j in new_gr]
            new_outs = [outs[j] for j in new_gr]
            new_time = [time[j] for j in new_gr]

            new_inds = inds.copy()
            new_inds.discard(prom_subgroups[i].cov_ind)

            print 'call sub_ident on split: %d' % (i)

            res_groups.extend(sub_ident(new_s, new_inds, covs, new_treat, new_outs, new_time,
                                        tree_h, recurs_level + 1, prom_subgroups[i]))

    print 'end call sub_ident at recurs_level: %d' % (recurs_level)

    return res_groups
"""

"""
def select_best(groups):
    best = []
    # 1 - left, 2 - right
    for gr in groups:
        if gr[1].get_p_val() < gr[2].get_p_val():
            best.append(tuple([gr[1], 1]))
        elif gr[1].get_p_val() > gr[2].get_p_val():
            best.append(tuple([gr[2], 2]))
        else:
            print 'Logrank p-values equal'
            raise
    return best
"""

"""
def select_best_subgroup(s, ind, covs, treat, outs, time):

    def get_cov_values(s, ind):
        return tuple(set([s[i][ind] for i in xrange(len(s))]))


#    def get_num_splits(s, ind, cov_type, values):
#        if cov_type == 'num':
#            return len(values)
#        else:
#            return 2 ** len(values) - 2


    def make_split(smp, ind, cov_type, values):
        if cov_type == 'num':
            for level in values:
                left, right = [], []
                for i in xrange(len(smp)):
                    if smp[i][ind] < level:
                        left.append(i)
                    else:
                        right.append(i)
                yield left, right, '>=%f' % (level), '<%f' % (level)
        else:
            for k in xrange(1, 2 ** len(values) - 1):
                left, right = [], []
                s = bin(k)[2:]
                if len(s) < len(values):
                    s = "".join(['0' for i in xrange(len(values) - len(s))]) + s
                level = tuple([values[i] for i in xrange(len(s)) if s[i] == '1'])
                not_level = tuple([values[i] for i in xrange(len(s)) if s[i] == '0'])
                for i in xrange(len(smp)):
                    if smp[i][ind] not in level:
                        left.append(i)
                    else:
                        right.append(i)
                yield left, right, 'in %s' % (str(level)), 'in %s' % (str(not_level))


    cov_type = covs[ind][1]
    cov_values = get_cov_values(s, ind)
    #n_split = get_num_splits(s, ind, covs, cov_values)
    treat_types = sorted(list(set(treat)))

    split_subgroup, split_z_val, split_p_logrank, split_pwr, split_better, split_levels = [], [], [], [], [], []
    split_p_crit = []
    split_surv_curves, split_n_surv, split_n_total = [], [], []

    for left, right, level, n_level in make_split(s, ind, cov_type, cov_values):
        if len(left) < dsc.minSubSize or len(right) < dsc.minSubSize:
            continue
        treat1, treat2 = [treat[i] for i in left], [treat[i] for i in right]
        outs1, outs2 = [outs[i] for i in left], [outs[i] for i in right]
        time1, time2 = [time[i] for i in left], [time[i] for i in right]

        t11, t12 = [i for i in xrange(len(treat1)) if treat1[i] == treat_types[0]], \
                   [i for i in xrange(len(treat1)) if treat1[i] == treat_types[1]]
        t21, t22 = [i for i in xrange(len(treat2)) if treat2[i] == treat_types[0]], \
                   [i for i in xrange(len(treat2)) if treat2[i] == treat_types[1]]

        s_curve11, n11_surv, n11_total = statTests.kaplan_meier([outs1[i] for i in t11], [time1[i] for i in t11])
        s_curve12, n12_surv, n12_total = statTests.kaplan_meier([outs1[i] for i in t12], [time1[i] for i in t12])
        s_curve21, n21_surv, n21_total = statTests.kaplan_meier([outs2[i] for i in t21], [time2[i] for i in t21])
        s_curve22, n22_surv, n22_total = statTests.kaplan_meier([outs2[i] for i in t22], [time2[i] for i in t22])

        z1, p_val1, pwr1, bet1 = statTests.logrank(treat1, outs1, time1, [s_curve11[-1], s_curve12[-1]])
        if z1 is None and p_val1 is None and pwr1 is None and bet1 is None:
            continue

        z2, p_val2, pwr2, bet2 = statTests.logrank(treat2, outs2, time2, [s_curve21[-1], s_curve22[-1]])
        if z2 is None and p_val2 is None and pwr2 is None and bet2 is None:
            continue

        if p_val1 > p_val2:
            best_group = left
            best_z = z1
            best_p_val = p_val1
            best_pwr = pwr1
            best_bet = bet1
            best_curves = (s_curve11, s_curve12)
            best_survivors = (n11_surv, n12_surv)
            best_total = (n11_total, n12_total)
        else:
            best_group = right
            best_z = z2
            best_p_val = p_val2
            best_pwr = pwr2
            best_bet = bet2
            best_curves = (s_curve21, s_curve22)
            best_survivors = (n21_surv, n22_surv)
            best_total = (n21_total, n22_total)

        split_subgroup.append(best_group)
        split_z_val.append(best_z)
        split_p_logrank.append(best_p_val)
        split_pwr.append(best_pwr)
        split_better.append(best_bet)
        split_surv_curves.append(best_curves)
        split_levels.append(level)
        split_n_surv.append(best_survivors)
        split_n_total.append(best_total)
        split_p_crit.append(split_criterion(z1,z2))

    left, right = select_num(s, ind, treat, outs, time) if cov_type == 'num' \
        else select_cat(s, ind, treat, outs, time)

    if len(left) < dsc.minSubSize or len(right) < dsc.minSubSize:
        pass

    treat1, treat2 = [treat[i] for i in left], [treat[i] for i in right]
    outs1, outs2 = [outs[i] for i in left], [outs[i] for i in right]
    time1, time2 = [time[i] for i in left], [time[i] for i in right]

    z1, p_val1, pwr1 = statTests.logrank(treat1, outs1, time1)
    z2, p_val2, pwr2 = statTests.logrank(treat2, outs2, time2)

    p = []
    for i in range(len(tests_left)):
        p.append(tuple([split_criterion(tests_left[i][0], tests_right[i][0]), i]))
    p.sort(key=lambda x: x[0], reverse=True)

    subgroups = []
    for i in range(len(p)):
        subgroups.append(spCl.BinarySplit().set_split(spCl.SplitContent().set_content()))
    wf.write_logrank()

    return split_subgroup, split_z_val, split_p_logrank, split_pwr, split_better, \
           split_surv_curves, split_n_surv, split_n_total, split_p_crit, split_levels
"""

"""
def get_values(s, inds):
    val_list = tuple([{} for i in inds])
    for i in inds:
        for k in range(len(s)):
            val_list[i].add(s[k][i])
    return val_list
"""


"""
def split_num(s, ind, treat, outs, time):

    def get_split_points(s, cov):
        vals = []
        for i in range(len(s)):
            vals.append(s[i][cov])
        else:
            return tuple(set(vals))

    levels = get_split_points(s, ind)

    for level in levels:
        left = right = []
        for i in range(len(s)):
            if s[i][ind] < level:
                left.append(i)
            else:
                right.append(i)

    return left, right
"""

"""
def select_cat(s, ind, treat, outs, time):

    def get_split_points(cov_values):
        points = []
        for k in range(1, 2 ** len(cov_values) - 1):
            s = bin(k)[2:]
            temp = []
            for i in range(len(s)):
                if s[i] == '1':
                    temp.append(set(cov_values[i]))
            points.append(temp)
        return tuple(points)

    levels = get_split_points(cov_values)
    n_splits = len(levels)

    tests_left, tests_right = []
    for k in range(n_splits):
        left = right = []
        for i in range(len(sample)):
            if sample[i][cov_index] not in levels[k]:
                left.append(sample[i])
            else:
                right.append(sample[i])
        if len(left) < dsc.minSubSize or len(right) < dsc.minSubSize:
            continue

        tests_left.append(tuple(statTests.logrank(left, titles, outs, tm, treat1, treat2) + levels[k]))
        tests_right.append(tuple(statTests.logrank(right, titles, outs, tm, treat1, treat2) + levels[k]))
            # T_left is a list of values of logrank test statistics for pre-specified treat_types
            # senslog_left is a list of values of sensitivities for each logrank test
    return tests_left, tests_right
"""

"""
def continuation_criterion(subgroups, par_sub, gamma=0.7):
    if par_sub is None:
        return range(len(subgroups))
    else:
        return [i for i in xrange(len(subgroups)) if subgroups[i].p_val_logrank <= gamma * par_sub.p_val_logrank]


def selection_criterion(x):
    pass
"""

"""
def sub_ident(sample, titles, cov_indices, outcomes, res_time, treatment):

    num_level = 1

    def insert_into_list(dest_list, src_list):
        return dest_list.extend(src_list).sort(key=spCl.sort_by_pvalue, reverse=True)[:dsc.numBestCov]

    # splitList is a list of result splits, after select_subgroup and, optionally, selection_criterion
    # cur_splitList is a list of splits given on

    n = len(list(set(treatment)))
    if n == 2:
        split_list = []

        if num_level == 1:
            for cov in cov_indices:
                val_list = get_values(sample, cov)
                if len(split_list):
                    split_list = insert_into_list(split_list, select_subgroup(sample, titles, cov, val_list))
                else:
                    split_list = select_subgroup(sample, titles, cov, val_list)

            split_list_div = continuation_criterion(split_list)
            if len(split_list_div):
                Q = queue.LifoQueue()
                k = len(split_list_div)
                for x in split_list_div:
                    Q.put(x)
                while Q.qsize():
                    cur_split = Q.get()
                    s = [sample[i] for i in cur_split.split_content.first] +\
                        [sample[i] for i in cur_split.split_content.second]





    else:
        split_list = [[] for x in range(int(n * (n-1) / 2))]

    while num_level < dsc.treeHeight:
        if num_level == 1:
            split_list = []
            for cov in cov_indices:
                    # promising subgroups
                    # extracting subgroups based on logrank test statistics and certain splitting criterion
                    val_list = get_values(sample, cov)
                    sp_l = select_subgroup(sample, titles, cov, val_list)
                    if len(split_list):

                    else:
                        split_list = sp_l

                    # whether we should split parent group into 2 child nodes or not
                    #child1 = continuation_criterion(child1)
                    #child2 = continuation_criterion(child2)
                    # candidate subgroups
                    # confirm statistically the significance of selected subgroups
                    #
                    subgroups = selection_criterion(subgroups)

        #numLevel += 1

    # apply subgroup selection procedure to splitList
    return split_list,
    #return subgroups
"""
