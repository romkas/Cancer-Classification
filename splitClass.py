from statTests import get_kmf_survival


class Subgroup(object):
    def __init__(self):
        self.group = None
        self.kmf1, self.kmf2 = None, None
        self.tt1, self.tt2 = None, None
        self.logrank = None
        self.pwr = None

        #self.sz = None
        #self.surv_n = None
        #self.surv = None
        #self.better = None

    def __setsz__(self):
        return len(self.group)

    def set_subgroup(self, gr, kmf1, tt1, kmf2, tt2, logrank, pwr):
        self.group = gr
        self.kmf1, self.kmf2 = kmf1, kmf2
        self.tt1, self.tt2 = tt1, tt2
        self.logrank = logrank
        self.pwr = pwr
        #self.sz = self.__setsz__()
        return

    def get_sz(self):
        return len(self.kmf1.durations), len(self.kmf2.durations)

    def get_surv(self):
        return get_kmf_survival(self.kmf1), get_kmf_survival(self.kmf2)

    def get_time_last_obs(self):
        return self.kmf1.timeline[-1], self.kmf2.timeline[-1]

"""
class Subgroup(object):
    def __init__(self):
        self.group = None
        self.cov_ind = None
        self.cov_name = None
        self.cov_level = None
        self.z = None
        self.p_val_logrank = None
        self.pwr = None
        self.p_val_crit = None
        self.sz = None
        self.surv_n = None
        self.surv = None
        self.better = None

        self.desc = None

    def set_subgroup(self, group, cov_ind, cov_name, cov_level, z, p_val_logrank, pwr,
                     p_val_crit, sz, surv_n, surv, better):
        self.group = group
        self.cov_ind = cov_ind
        self.cov_name = cov_name
        self.z = z
        self.p_val_logrank = p_val_logrank
        self.pwr = pwr
        self.p_val_crit = p_val_crit
        self.sz = sz
        self.surv_n = surv_n
        self.surv = surv
        self.better = better

    def set_descr(self, desc):
        self.desc = desc
"""

"""
class BinarySplit(object):
    def __init__(self):
        #self.patients = None
        self.pats_first = None
        self.pats_second = None
        self.cov_index = None
        self.cov_level = None
        #self.cov_name = None
        self.lrank_val = None
        self.p_val_lrank = None
        self.power = None
        self.p_val_split = None
        self.t_names = None
        self.t_better = None
    
    def __len__(self):
        return len(self.patients)

    def __repr__(self):
        s = 'Size : ( %d, %d )\nTreat : %s, %s\nBetter : %d\nCovar : %s' %\
            (len(self.pats_first), len(self.pats_second),
             self.t_names[0], self.t_names[1], self.t_better, self.cov_name)
        if Dsc.covTypes[self.cov_index] == 'num':
            s += 'Levels : %d' % self.cov_level
        else:
            s += 'Levels : '
            for x in self.cov_level:
                s += str(x) + ', '
        s += '\nLogrank : %f\nPower : %f\nP-value : %f\n' %\
             (self.p_val_lrank, self.power, self.p_val_split)
        return s

    def __lt__(self, other):
        return self.p_val_split < other.p_val_split

    def __le__(self, other):
        return self.p_val_split <= other.p_val_split

    def __gt__(self, other):
        return self.p_val_split > other.p_val_split

    def __ge__(self, other):
        return self.p_val_split >= other.p_val_split

    def __eq__(self, other):
        return self.p_val_split == other.p_val_split

    def __ne__(self, other):
        return self.p_val_split != other.p_val_split
    
    def set_content(self, patients=None, pats_first=None, pats_second=None, cov_index=-1, cov_level=-1, cov_name='',
                    lrank_val=-1, p_val_lrank=-1, power=-1,
                    p_val_split=-1, t_names=('t1', 't2'), t_better=-1):
        self.patients = patients
        self.pats_first = pats_first
        self.pats_second = pats_second
        self.cov_index = cov_index
        self.cov_level = cov_level
        self.cov_name = cov_name
        self.lrank_val = lrank_val
        self.p_val_lrank = p_val_lrank
        self.power = power
        self.p_val_split = p_val_split
        self.t_names = t_names
        self.t_better = t_better
        return self

    def set_split(self, patients, pats_first, pats_second, cov_index, cov_level, cov_name):
        self.patients, self.pats_first, self.pats_second = pats_first, pats_second, patients
        self.cov_index, self.cov_level, self.cov_name = cov_index, cov_level, cov_name
        return self

    def set_logrank(self, lrank_val, p_val, power):
        self.lrank_val = lrank_val
        self.p_val_lrank = p_val
        self.power = power
        return self

    def set_p_val_split(self, p):
        self.p_val_split = p
        return self

    def set_better(self, t_names):
        self.t_names = t_names
        if self.p_val_lrank > 0:
            self.t_better = 1
        elif self.p_val_lrank == 0:
            self.t_better = 0
        else:
            self.t_better = -1
        return self

    def set_content_ws(self, spl):
        self.set_split(spl.patients, spl.pats_first, spl.pats_second, spl.cov_index, spl.cov_level, spl.cov_name)
        self.set_logrank(spl.lrank_val, spl.p_val_lrank, spl.power)
        self.set_p_val_split(spl.p_val_split)
        self.set_better(spl.t_names)
        return self

    def get_patients(self):
        return self.patients
"""

"""
class Node(object):
    def __init__(self):
        self.split_content = None  # content of BinarySplit class
        self.parent = None  # reference to parent split
        self.left = None    # left child
        self.right = None   # right child

    def set_parent(self, p_par):
        self.parent = p_par
        return self

    def set_left(self, p_left):
        self.left = p_left
        return self

    def set_right(self, p_right):
        self.right = p_right
        return self

    def set_split(self, data, p_par, p_left, p_right):
        self.split_content.set_content_ws(data)
        self.parent = p_par
        self.left = p_left
        self.right = p_right
        return self

    def get_p_val(self):
        return self.split_content.p_val_lrank

    def get_group(self, s):
        return [s[i] for i in self.split_content.get_patients()]


def sort_by_p_val(x):
    return x.p_val_crit
"""