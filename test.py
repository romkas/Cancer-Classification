from subgroupIdentify import subgroup_identification as subid
from plotting import plot_bin_tree
from main import select_binary_groups
from readfile import load_data2002


fname1 = 'D:\\Users\\ROM-VAIO\\Copy\\CoursePaper3\\DATA\\all-2002-rg1-6-60.txt'
fname2 = 'D:\\Users\\ROM-VAIO\\Copy\\CoursePaper3\\DATA\\all-2002-rg2-6-60.txt'
a_logrank = 0.01
cov_at_level = 3

s, nobs, ncov = load_data2002(fname2)

g = subid(s,a_logrank,cov_at_level)

trees = select_binary_groups(g)

for t in trees:
    plot_bin_tree(t, node_sz=600)
