from matplotlib import pyplot as plt
from Queue import Queue as Q
from numpy import arange
from networkx import draw_networkx, draw_networkx_edge_labels

def plot_kmf(group, dsc):
    def get_pos_textbox(group):
        max_obs = max(group.get_time_last_obs())
        min_surv = min(group.get_surv())
        x = 0.7 * max_obs
        y = 1 - 0.35 * (1 - min_surv)
        return x, y

    plt.figure()
    ax = plt.subplot(111)
    group.kmf1.plot(ax=ax, ci_show=False)
    group.kmf2.plot(ax=ax, ci_show=False)
    ax.set_title(dsc, bbox={'facecolor': 'white', 'pad': 5})
    s = 'p-value: %f\npower:  %f' % (group.res_logrank.p_value, group.pwr)
    x, y = get_pos_textbox(group)
    ylim_top, ylim_bot = 1, max(0, min(group.get_surv()) - 0.1)
    ax.text(x, y, s, bbox={'facecolor': 'white', 'pad': 7})
    ax.set_ylim(ylim_bot, ylim_top)


def make_layout(groups, x_range, y_range, root='0'):
    def get_graph_params(groups):
        num_levels = 1
        max_nodes_at_level = 1
        level = Q()
        level.put(root)
        while True:
            temp = []
            while level.qsize() > 0:
                temp.extend(groups.successors(level.get()))
            if temp:
                num_levels += 1
                if max_nodes_at_level < len(temp):
                    max_nodes_at_level = len(temp)
                while not temp:
                    level.put(temp.pop())
            else:
                break
        return num_levels, max_nodes_at_level

    def get_x_tick():
        pass

    def get_y_tick(y_range, nlevels):
        step = (y_range[1] - y_range[0]) / (nlevels + 2)
        return arange(y_range[0], y_range[1], step)[1:]

    x_span, y_span = float(x_range[1] - x_range[0]), float(y_range[1] - y_range[0])
    nlevels, natlvl = get_graph_params(groups)
    x_tick, y_tick = get_x_tick(), get_y_tick(y_range, nlevels)
    pos = {root: (float(x_span) / 2, float(y_span))}
    for nlvl in ''.join([str(x) for x in xrange(1, nlevels)]):



def plot_tree(groups):
    plt.figure()
    ax = plt.subplot(111)
    x_range, y_range = (0, 10), (0, 10)
    root = '0'
    pos = make_layout(groups, x_range, y_range, root=root)
