from readfile import load_data2002
from main import select_candidates as secan
from main import permute_sample
from subgroupIdentify import subgroup_identification as subid


fname2 = 'D:\\Copy\\CoursePaper3\\DATA\\all-2002-rg2-6-60.txt'
a0 = 0.05
s2, nobs2, ncov2 = load_data2002(fname2)
gr_file2_p = subid(permute_sample(s2), a0)
cands_file2 = secan(gr_file2_p)
for c in cands_file2:
    print c
    print gr_file2_p.node[c]['content'].pwr
    print gr_file2_p.node[c]['content'].res_logrank
