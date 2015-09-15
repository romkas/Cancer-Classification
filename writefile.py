def write_subgroups(smp, covs, groups, desc_t, desc_c, surv, desc_common):
    f1, f2, f3, f4 = open(desc_t, 'a'), open(desc_c, 'a'), open(surv, 'a'), open(desc_common, 'a')
    f2.write('Size;Sz1;Sz2;Num_surv1;Num_surv2;Covar;Level;Bet_tr;Z;P-val;Power' + '\n')
    for gr in groups:
        f1.write(gr.desc + '\n\n')

        f2.write('%d;%d;%d;%d;%d;%s;%s;%d;%.4f;%.4f;' % (gr.sz[0] + gr.sz[1], gr.sz[0], gr.sz[1], gr.surv_n[0],
                                                         gr.surv_n[1], gr.cov_name, str(gr.cov_level), gr.better,
                                                         gr.z, gr.p_val_logrank))
        if gr.pwr is None:
            f2.write('None\n')
        else:
            f2.write('%.4f\n' % (gr.pwr))

        f3.write(';'.join([str(x) for x in gr.surv[0]]) + '\n' + ';'.join([str(x) for x in gr.surv[1]]) + '\n\n')

        pos_treat = covs.index(('Rand1', 'null'))
        pos_tod = covs.index(('Tod', 'null'))
        pos_time = covs.index(('Time', 'null'))

        f4.write('Treat;Tod;Time\n')
        for i in gr.group:
            s = ';'.join([str(smp[i][pos_treat]), str(smp[i][pos_tod]), str(smp[i][pos_time])])
            f4.write(s + '\n')
        f4.write('\n')

    f1.close()
    f2.close()
    f3.close()
    f4.close()

    return


"""
def write_logrank(T, pvalue, sens):
    fname = 'D:\\Copy\\CoursePaper3\\out-files\\logrank.txt'
    with open(fname, 'w+') as outfile:
        for i in range(len(T)):
            outfile.write('\t'.join([str(T[i]), str(pvalue[i]), str(sens[i])]))
            outfile.write('\n\n')

    return 0
"""
