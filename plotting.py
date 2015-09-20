from matplotlib import pyplot as plt
from Queue import Queue
from networkx import draw_networkx, draw_networkx_edge_labels
import desc as dsc

def plot_kmf(group, desc):
    def get_pos_textbox(grp):
        max_obs = max(grp.get_time_last_obs())
        min_surv = min(grp.get_surv())
        x = 0.7 * max_obs
        y = 1 - 0.35 * (1 - min_surv)
        return x, y

    plt.figure()
    ax = plt.subplot(111)
    group.kmf1.plot(ax=ax, ci_show=False)
    group.kmf2.plot(ax=ax, ci_show=False)
    ax.set_title(desc, bbox={'facecolor': 'white', 'pad': 5})
    s = 'p-value: %f\npower:  %f' % (group.logrank.p_value, group.pwr)
    x_coord, y_coord = get_pos_textbox(group)
    ylim_top, ylim_bot = 1, max(0, min(group.get_surv()) - 0.1)
    ax.text(x_coord, y_coord, s, bbox={'facecolor': 'white', 'pad': 7})
    ax.set_ylim(ylim_bot, ylim_top)


def make_layout(groups, node_size, root='0'):
    def calc_pos_current_node(pos_parent, node_sz, sd):
        # sd == -1 if left; == 1 if right
        def calc_step_down(nd_sz):
            return 2

        def calc_step_side(nd_sz):
            return 1

        step_down = calc_step_down(node_sz)
        step_side = calc_step_side(node_sz)
        x_coord = pos_parent[0] - step_side if sd < 0 else pos_parent[0] + step_side
        return x_coord, pos_parent[1] - step_down

    def node_location(nd):
        return -1 if nd.split('-')[2] == 'left' else 1

    pos = {root: (0, 0)}
    q = Queue()
    q.put(root)
    while not q.empty():
        node = q.get()
        scs = groups.successors(node)
        for nd in scs:
            side = -1 if node_location(nd) < 0 else 1
            pos[nd] = calc_pos_current_node(pos[node], node_sz=node_size, sd=side)
            q.put(nd)
    return pos


def get_edge_label(groups, u, v):
    return groups.get_edge_data(u, v)[dsc.edge_content]


def make_edge_labels(groups):
    e = groups.edges()
    values = [get_edge_label(groups, u, v) for u, v in e]
    labels = [(e[i], values[i]) for i in xrange(len(e))]
    return {key: val for (key, val) in labels}


def plot_bin_tree(groups, node_size):
    plt.figure()
    ax = plt.subplot(111)
    root = '0'
    pos = make_layout(groups, node_size=node_size, root=root)
    edge_labels = make_edge_labels(groups)
    draw_networkx(groups, pos=pos, ax=ax, with_labels=True, node_size=node_size, node_shape='s')
    draw_networkx_edge_labels(groups, pos=pos, ax=ax, edge_labels=edge_labels)
