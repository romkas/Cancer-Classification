import desc


def load_data(fname):
    f = open(fname, 'r')
    sample = [x.replace(',', '.').split('\t') for x in f.readlines()]
    f.close()
    return sample


def filter_data(sample, min_time_research=desc.min_time_research):
    s = []
    s_dict = {}
    ind_time = 2
    for i in xrange(len(sample)):
        a = [float(x) for x in sample[i]]
        if a[ind_time] >= min_time_research:
            s.append(a)
    for i, cov_name in enumerate(desc.covariates):
        s_dict[cov_name] = [s[k][i] for k in xrange(len(s))]
    return s_dict
