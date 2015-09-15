import Descriptions as dsc


def load_data2002(fname):
    f = open(fname, 'r')
    sample = [x.replace(',', '.').split('\t') for x in f.readlines()]
    f.close()
    if not f.closed:
        exit('file is not closed')
    s = filter_data(sample)
    return s, len(s['Rand']), len(s) - 3


def filter_data(sample):
    s = []
    s_dict = {}
    for i in xrange(len(sample)):
        a = [float(x) for x in sample[i]]
        if a[dsc.covariates2002.index('Time')] >= dsc.min_time_research:
            s.append(a)
    for i in enumerate(dsc.covariates2002):
        s_dict[i[1]] = [s[k][i[0]] for k in xrange(len(s))]
    return s_dict


def dict_convert(sample):
    s = {}
    if type(sample) == list:
        for item in enumerate(dsc.covariates2002):
            s[item[1]] = [sample[k][item[0]] for k in xrange(len(sample))]
        sample = dict(s)
    else:
        exit('to_dict() type error')
    return sample


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


def remove_min_time(sample):
    if type(sample) == dict:
        n = len(sample[dsc.covariates2002[0]])
    else:
        exit('remove_min_time() error')
    inds = [i for i in xrange(n) if sample['Time'][i] >= dsc.min_time_research]
    for cov in sample:
        sample[cov] = [sample[cov][k] for k in inds]
    return sample