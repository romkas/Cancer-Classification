import desc as dsc


def load_data2002(fname):
    f = open(fname, 'r')
    sample = [x.replace(',', '.').split('\t') for x in f.readlines()]
    f.close()
    if not f.closed:
        exit('file is not closed')
    return sample


def filter_data(sample, min_time_research=dsc.min_time_research):
    s = []
    s_dict = {}
    ind_time = 2
    for i in xrange(len(sample)):
        a = [float(x) for x in sample[i]]
        if a[ind_time] >= min_time_research:
            s.append(a)
    for i, cov_name in enumerate(dsc.covariates2002):
        s_dict[cov_name] = [s[k][i] for k in xrange(len(s))]
    return s_dict

"""
def dict_convert(sample):
    s = {}
    if type(sample) == list:
        for i, cov_name in enumerate(dsc.covariates2002):
            s[cov_name] = [sample[k][i] for k in xrange(len(sample))]
        sample = dict(s)
    else:
        exit('to_dict() type error')
    return sample
"""
"""
def to_float(sample):
    if type(sample) == dict:
        for cov in sample.keys():
            sample[cov] = [float(x) for x in sample[cov]]
    elif type(sample) == list:
        for i in xrange(len(sample)):
            sample[i] = [float(x) for x in sample[i]]
    else:
        exit('to_float() type error')
    return sample
"""
"""
def remove_min_time(sample):
    rand = dsc.covariates2002[0]
    time = dsc.covariates2002[2]
    n = len(sample[rand])
    inds = [i for i in xrange(n) if sample[time][i] >= dsc.min_time_research]
    for cov in sample:
        sample[cov] = [sample[cov][k] for k in inds]
    return sample
"""