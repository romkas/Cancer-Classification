import desc
import splitClass as sc
import statTests as stat
from math import fabs
from math import sqrt
from scipy.stats import norm
from networkx import DiGraph
from numpy.random import permutation as pm
from numpy import arange


global_counter = 0


def subgroup_identification(sample, mode,
                            a_logrank=desc.a_logrank,
                            cov_at_level=desc.cov_at_level,
                            min_sub_size=desc.min_sub_size):
    subgroups = DiGraph()
    cov_used = {key: val for key in sample.keys() for val in [False]}
    subgroups.add_node(desc.root)
    global global_counter
    global_counter = 1
    sub_ident(sample, subgroups, cov_used, a_logrank, cov_at_level, min_sub_size, mode,
              parent=desc.root, recurs_level=1)
    return subgroups


def get_group(groups, node_name):
    return groups.node[node_name][desc.node_content]


def sub_ident(sample, subgroups, cov_used, a_logrank, cov_at_level, min_sub_size, mode, parent, recurs_level):
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

    def create_node_name(rec_lvl, mod):
        if not (mod == 'left' or mod == 'right'):
            exit('create_node_name() error')
        return '-'.join([str(rec_lvl), str(global_counter), mod])

    def better_group(grp1, grp2):
        return 1 if grp1.logrank.p_value < grp2.logrank.p_value else -1

    def select_subsample(samp, group_inds):
        subsamp = {key: val for key in samp.keys() for val in [[]]}
        for cov in subsamp.keys():
            subsamp[cov] = [samp[cov][i] for i in xrange(len(samp[cov])) if i in group_inds]
        return subsamp

    group0 = []
    rand, tod, time = desc.covariates[0], desc.covariates[1], desc.covariates[2]
    for cov_name, val in cov_used.items():
        if val or cov_name == rand or cov_name == tod or cov_name == time:
            continue
        group = select_best_k_splits(sample, cov_name, a_logrank, cov_at_level, min_sub_size, mode)
        merge(group, group0)
    gamma = select_cont_param()
    better = []
    groups_next = []
    cov_next = []
    global global_counter
    for split in group0:
        name1 = create_node_name(recurs_level, 'left')
        global_counter += 1
        name2 = create_node_name(recurs_level, 'right')
        global_counter += 1
        subgroups.add_node(name1, {desc.node_content: split[0]})
        subgroups.add_edge(parent, name1, {desc.edge_content: split[1]})
        subgroups.add_node(name2, {desc.node_content: split[2]})
        subgroups.add_edge(parent, name2, {desc.edge_content: split[3]})
        better.append(
            name1 if better_group(get_group(subgroups, name1), get_group(subgroups, name2)) > 0 else name2
        )
        if parent == desc.root:
            groups_next.append(better[-1])
            cov_next.append(split[4])
        else:
            if continuation_criterion(get_group(subgroups, better[-1]), get_group(subgroups, parent), gamma):
                groups_next.append(better[-1])
                cov_next.append(split[4])
    for i in xrange(len(groups_next)):
        subsample = select_subsample(sample, get_group(subgroups, groups_next[i]).group)
        new_cov_used = cov_used.copy()
        new_cov_used[cov_next[i]] = True
        sub_ident(subsample, subgroups, new_cov_used, a_logrank, cov_at_level, min_sub_size, mode,
                  groups_next[i], recurs_level + 1)
    return subgroups


def continuation_criterion(current, parent, gamma):
    return current.logrank.p_value <= gamma * parent.logrank.p_value


def select_cont_param():
    return 1


def select_best_k_splits(sample, cov_name, a_logrank, cov_at_level, min_sub_size, mode):
    # returning values:
    # group = [ tuple(group1, desc1, group2, desc2, cov_name, p_value_split), ... ]
    def get_cov_values(samp, cov):
        return tuple(set(samp[cov]))

    def make_split(samp, cov, vals, mod):
        if cov == 'Sex':
            l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                   [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
            yield l, r, 'Male', 'Female'
        elif cov == 'Immun':
            if mod == desc.modes[0]:
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
                    yield l, r, 'Immun-B in %s' % (str(level)), 'Immun-B in %s' % (str(not_level))
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
                    yield l, r, 'Immun-T in %s' % (str(level)), 'Immun-T in %s' % (str(not_level))
            elif mod == desc.modes[1] or mode == desc.modes[2]:
                set_vals = set(vals)
                level = [{1, 2, 3, 4}, {10}, {12}, {13}, {5, 6, 7, 8}, {11}, {14, 15, 16, 17, 18, 19, 20}]
                not_level = [set_vals.difference(lvl) for lvl in level]
                for k in xrange(len(level)):
                    l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] in not_level[k]],\
                           [j for j in xrange(len(samp[cov])) if samp[cov][j] in level[k]]
                    yield l, r, 'Immun in %s' % ('{' + str(not_level[k])[5:-2] + '}'),\
                                'Immun-B in %s' % ('{' + str(level[k])[5:-2] + '}') if k < 4\
                                else 'Immun-T in %s' % ('{' + str(level[k])[5:-2] + '}')
            else:
                exit('make_split() error')
        elif cov == 'CNS':
            if mod == desc.modes[0]:
                l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                       [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
                yield l, r, 'CNS = 1', 'CNS = 2'
            elif mod == desc.modes[1] or mod == desc.modes[2]:
                for lvl in vals:
                    l, r = [], []
                    for j in xrange(len(samp[cov])):
                        if samp[cov][j] < lvl:
                            l.append(j)
                        else:
                            r.append(j)
                    level = [x for x in vals if x < lvl]
                    not_level = [x for x in vals if x >= lvl]
                    yield l, r, '%s = %s' % (cov, '{' + str(level)[1:-1] + '}'),\
                                '%s = %s' % (cov, '{' + str(not_level)[1:-1] + '}')
            else:
                exit('make_split() error')
        elif cov == 'Mediastinum':
            if mod == desc.modes[0] or mod == desc.modes[1]:
                l, r = [j for j in xrange(len(samp[cov])) if samp[cov][j] == 1],\
                       [j for j in xrange(len(samp[cov])) if samp[cov][j] == 2]
                yield l, r, 'Mediastinum = 1', 'Mediastinum = 2'
            elif mod == desc.modes[2]:
                for lvl in vals:
                    l, r = [], []
                    for j in xrange(len(samp[cov])):
                        if samp[cov][j] < lvl:
                            l.append(j)
                        else:
                            r.append(j)
                    level = [x for x in vals if x < lvl]
                    not_level = [x for x in vals if x >= lvl]
                    yield l, r, '%s = %s' % (cov, '{' + str(level)[1:-1] + '}'),\
                                '%s = %s' % (cov, '{' + str(not_level)[1:-1] + '}')
            else:
                exit('make_split() error')
        elif cov == 'Age' or cov == 'Leuc' or cov == 'Leber' or cov == 'Milz':
            for level in vals:
                l, r = [], []
                for j in xrange(len(samp[cov])):
                    if samp[cov][j] < level:
                        l.append(j)
                    else:
                        r.append(j)
                yield l, r, '%s < %f' % (cov, level), '%s >= %f' % (cov, level)
        else:
            exit('make_split() error')

    def convert_data(outcm):
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

    treatment = sample[desc.covariates[0]]  # Rand - protocol
    outcome = sample[desc.covariates[1]]  # Tod - alive, lost, dead
    time = sample[desc.covariates[2]]  # Time - lifetime
    treatment_type = sorted(list(set(treatment)))
    best_k_splits = []
    cov_values = get_cov_values(sample, cov_name)
    if len(cov_values) == 1:
        return best_k_splits
    for left, right, desc_left, desc_right in make_split(sample, cov_name, cov_values, mode):
        if len(left) < min_sub_size or len(right) < min_sub_size:
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


def permute_sample(sample):
    rand = desc.covariates[0]
    tod = desc.covariates[1]
    time = desc.covariates[2]
    nobs = len(sample[rand])
    inds = pm(nobs)
    s = {key: range(nobs) for key in sample.keys()}
    s[rand] = list(sample[rand])
    s[tod] = list(sample[tod])
    s[time] = list(sample[time])
    for i in xrange(len(inds)):
        for cov in s.keys():
            if cov == rand or cov == tod or cov == time:
                continue
            s[cov][i] = sample[cov][inds[i]]
    return s


def select_candidates(groups):
    return [node for node in groups.nodes_iter() if not groups.successors(node)]


def print_group(groups, node_name, f):
    f.write(node_name + '\n')
    n1, n2 = get_group(groups, node_name).get_sz()
    f.write('Left: %d\tRight: %d\n' % (n1, n2))
    f.write('Power: %f\n' % get_group(groups, node_name).pwr)
    f.write(str(get_group(groups, node_name).logrank) + '\n')


def resampling(sample,
               upper, lower, step,
               alpha0, cov_at_level, min_sub_size=desc.min_sub_size,
               out_file=None,
               n=100):
    cutoffs = arange(upper, lower, step)
    max_cutoff = 1
    flag = 1 if out_file is not None else 0
    for cut in cutoffs:
        print 'Cutoff = %f' % cut
        count = 0
        print 'Run: ',
        for i in xrange(n):
            print '%d' % i,
            groups = subgroup_identification(permute_sample(sample), cut, cov_at_level, min_sub_size)
            candidates = select_candidates(groups)
            candidates.sort()
            if flag:
                f = open(out_file, 'a')
            for cs in candidates:
                if flag:
                    print_group(groups, cs, f)
                if get_group(groups, cs).logrank.is_significant:
                    count += 1
                    break
            if flag:
                f.write('\n')
                f.close()
            if count > n * alpha0:
                break
        print ''
        if count <= n * alpha0:
            max_cutoff = cut
            break
    return max_cutoff


def select_binary_groups(groups, root=desc.root):
    def get_bin_group(grps, bin_grps, rt, graph):
        level = grps.successors(rt)
        if not level:
            bin_grps.append(graph)
            return bin_grps
        elif len(level) == 2:
            graph.add_nodes_from(level)
            graph.add_edge(rt, level[0], grps.get_edge_data(rt, level[0]))
            graph.add_edge(rt, level[1], grps.get_edge_data(rt, level[1]))
            if grps.successors(level[0]):
                get_bin_group(grps, bin_grps, rt=level[0], graph=graph)
            else:
                get_bin_group(grps, bin_grps, rt=level[1], graph=graph)
        else:
            level.sort()
            for i in xrange(0, len(level), 2):
                g = DiGraph(graph)
                g.add_nodes_from([level[i], level[i + 1]])
                g.add_edge(rt, level[i], grps.get_edge_data(rt, level[i]))
                g.add_edge(rt, level[i + 1], grps.get_edge_data(rt, level[i + 1]))
                if grps.successors(level[i]):
                    get_bin_group(grps, bin_grps, rt=level[i], graph=g)
                else:
                    get_bin_group(grps, bin_grps, rt=level[i + 1], graph=g)

    dg = DiGraph()
    dg.add_node(root)
    bin_groups = []
    get_bin_group(groups, bin_grps=bin_groups, rt=root, graph=dg)
    return bin_groups
