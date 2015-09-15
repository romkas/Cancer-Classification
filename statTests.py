#import xalglib as xa
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats.distributions import norm
from math import log
from math import sqrt
"""
def logrank(sample, titles, treat_types):
    posTime = titles.index('Time')
    posTod = titles.index('Tod')
    posTreat = titles.index('Rand1')

    min_time = min([[posTime] for x in sample])
    max_time = max([[posTime] for x in sample])

    if len(treat_types) == 2:
        T1, T2 = treat_types[0], treat_types[1]

        s_treat1, s_treat2 = [x for x in sample if x[posTreat] == T1],\
                             [x for x in sample if x[posTreat] == T2]

        n1, n2 = len(s_treat1), len(s_treat2)

        U, Var_U = [0], [0]

        for t in range(min_time, max_time + 1):
            current1, current2 = [x[posTod] for x in s_treat1 if x[posTime] == t], \
                                [x[posTod] for x in s_treat2 if x[posTime] == t]

            n1_death, n2_death = len([x for x in current1 if x]), len([x for x in current2 if x])

            if not n1_death and not n2_death:
                continue

            n1_death_expec, n2_death_expec = n1 * (n1_death + n2_death) / (n1 + n2), \
                                            n2 * (n1_death + n2_death) / (n1 + n2)

            U[0] += n1_death - n1_death_expec
            Var_U[0] += n1*n2*(n1_death+n2_death)*(n1+n2-n1_death-n2_death)/(n1+n2)/(n1+n2)/(n1+n2-1)

            n1 -= n1_death
            n2 -= n2_death

    elif len(treat_types) > 2:
        pass
    else:
        pass

    sens = logrank_sens(treat_types) # a list of sensitivities for each logrank test

    return [U[i] * U[i] / Var_U[i] for i in range(len(U))], sens
"""

def kaplan_meier(out, t, ttype):
    def make_label(ttype, nobs):
        return 'Rand%d; %d obs.' % (ttype, nobs)

    kmf = KaplanMeierFitter()
    kmf.fit(t, event_observed=out, label=make_label(ttype=ttype, nobs=len(out)))
    return kmf


def logrank(out1, time1, out2, time2, alpha=0.05):
    return logrank_test(time1, time2, out1, out2, alpha=1-alpha)


def logrank_power(n, surv1, surv2, alpha=0.05):
        d = n * (2 - surv1 - surv2)
        if surv1 == 1 or surv2 == 1:
            return 0
        elif surv1 == 0 or surv2 == 0:
            return -1
        phi = log(surv1) / log(surv2) if surv1 < surv2 else log(surv2) / log(surv1)
        z_a = norm.ppf(1 - alpha)
        z_1_beta = sqrt(d * (1 - phi) * (1 - phi) / (1 + phi) / (1 + phi)) - z_a
        return norm.cdf(z_1_beta)

"""
def logrank(treat, outs, time, surv):

    def get_mode():
        return 1

    def get_outcome():
        return [2]

    def count_outcomes(outs):
        return len([x for x in outs if x in get_outcome()])

    def logrank_test(x, mode=get_mode()):
        if mode == 1:
            return xa.chisquaredistribution(1, x)
        else:
            return xa.normaldistribution(x)

    def logrank_power(n1, n2, surv1, surv2):
        alpha = 0.05

        n = min(n1, n2)
        d = n * (2 - surv1 - surv2)
        phi = math.log(surv2) / math.log(surv1)
        z_a = xa.invnormaldistribution(1 - alpha)
        z_1_beta = math.sqrt(d * (1 - phi) * (1 - phi) / (1 + phi) / (1 + phi)) - z_a

        return xa.normaldistribution(z_1_beta)

    def calc_z_logrank(u, var, mode=get_mode()):
        if mode == 1:
            return float(u) * u / var
        else:
            return float(u) / math.sqrt(var)

    min_time = min(time)
    max_time = max(time)

    types = tuple(set(treat))
    n1 = len([x for x in treat if x == types[0]])
    n2 = len(treat) - n1

    if surv[0] >= surv[1]:
        surv1, surv2, bet = surv[0], surv[1], 1
    else:
        surv1, surv2, bet = surv[1], surv[0], 2

    u, var_u = 0, 0
    n1_alive, n2_alive = n1, n2

    if surv1 == 1 and surv2 == 1 or surv1 == 0 and surv2 == 0:
        return None, None, None, None

    for t in xrange(int(min_time), int(max_time) + 1):
        group1 = [outs[i] for i in range(n1) if time[i] == t]
        group2 = [outs[i] for i in range(n2) if time[i] == t]

        n1_dead = count_outcomes(group1)
        n2_dead = count_outcomes(group2)

        if not n1_dead and not n2_dead:
            continue

        n1_dead_expec = float(n1_alive) * (n1_dead + n2_dead) / (n1_alive + n2_alive)
        # n2_death_expec = n2 * (n1_death + n2_death) / (n1 + n2)

        u += n1_dead - n1_dead_expec
        var_u += float(n1_alive) * n2_alive * (n1_dead + n2_dead) * (n1_alive + n2_alive - n1_dead - n2_dead) / \
                 (n1_alive + n2_alive) / (n1_alive + n2_alive) / (n1_alive + n2_alive-1)

        n1_alive -= n1_dead
        n2_alive -= n2_dead

    z = calc_z_logrank(u, var_u)
    p_val = logrank_test(z)

    pwr = logrank_power(n1, n2, surv1, surv2) if surv1 not in {0, 1} and surv2 not in {0, 1} else None

    return z, p_val, pwr, bet
"""

"""
def kaplan_meier(outs, time):

    def count_outcomes(outs, time, t, *args):
        return len([x for x in [outs[i] for i in xrange(len(time)) if time[i] == t] if x in args])

    min_time = min(time)
    max_time = max(time)
    surv = 1
    s_curv = []
    n_surv_cur = len(outs)
    for t in xrange(int(min_time), int(max_time) + 1):
        n_death = count_outcomes(outs, time, t, 2)
        surv *= (1 - float(n_death) / n_surv_cur)
        s_curv.append(surv)
        n_surv_cur -= n_death

    return s_curv, n_surv_cur, len(outs)
"""