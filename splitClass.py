from statTests import get_kmf_survival


class Subgroup(object):
    def __init__(self):
        self.group = None
        self.kmf1, self.kmf2 = None, None
        self.tt1, self.tt2 = None, None
        self.logrank = None
        self.pwr = None
        # self.sz = None
        # self.surv_n = None
        # self.surv = None
        # self.better = None

    def __setsz__(self):
        return len(self.group)

    def set_subgroup(self, gr, kmf1, tt1, kmf2, tt2, logrank, pwr):
        self.group = gr
        self.kmf1, self.kmf2 = kmf1, kmf2
        self.tt1, self.tt2 = tt1, tt2
        self.logrank = logrank
        self.pwr = pwr
        # self.sz = self.__setsz__()
        return

    def get_sz(self):
        return len(self.kmf1.durations), len(self.kmf2.durations)

    def get_surv(self):
        return get_kmf_survival(self.kmf1), get_kmf_survival(self.kmf2)

    def get_time_last_obs(self):
        return self.kmf1.timeline[-1], self.kmf2.timeline[-1]
