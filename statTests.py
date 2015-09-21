from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats.distributions import norm
from math import log
from math import sqrt


def get_kmf_survival(kmf):
    return kmf.survival_function_.get_values()[-2][0]


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
