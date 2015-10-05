from matplotlib import pyplot as plt
from Queue import Queue
from networkx import draw_networkx_nodes, draw_networkx_edges, draw_networkx_edge_labels
from math import pow
import desc as dsc


def get_edge_label(groups, u, v):
    return unicode(groups.get_edge_data(u, v)[dsc.edge_content])


def get_description(groups, node, root=dsc.root):
    s = u''
    count = 0
    f = True
    while node != root:
        parent = groups.predecessors(node)[0]
        if parent == root:
            s += get_edge_label(groups, parent, node)
        else:
            s += u''.join([get_edge_label(groups, parent, node), u';'])
        count += 1
        node = parent
        if f and count >= 3:
            if node != root:
                s += u'\n'
                f = False
    return s


def plot_kmf(group, name, desc):
    def get_pos_textbox(grp):
        # max_obs = max(grp.get_time_last_obs())
        # min_surv = min(grp.get_surv())
        # x = 0.7 * max_obs
        # y = 1 - 0.35 * (1 - min_surv)
        # return x, y
        return 0.01, 0.98  # relative coords

    plt.figure(name)
    ax = plt.subplot(111)
    group.kmf1.plot(ax=ax, ci_show=False, linestyle='--', color='k')
    group.kmf2.plot(ax=ax, ci_show=False, color='k')
    ax.set_title(desc, bbox={'facecolor': 'white', 'pad': 5})
    ax.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
    s = u'p-value: %f\npower:   %f' % (group.logrank.p_value, group.pwr)
    x_coord, y_coord = get_pos_textbox(group)
    max_timeline = max(group.get_time_last_obs())
    xlim_right, xlim_left = 1.01 * max_timeline, 0
    ylim_top, ylim_bot = 1 + 0.23, min(group.get_surv()) - 0.1
    ax.text(x_coord, y_coord, s, bbox={'facecolor': 'white', 'pad': 5},
            verticalalignment='top', horizontalalignment='left', transform=ax.transAxes)
    # ax.annotate(s, xy=(x_coord, y_coord), xycoords='axes fraction',
    #             xytext=(x_coord, y_coord), textcoords='axes fraction',
    #             bbox={'facecolor': 'white', 'pad': 7})
    ax.set_xlim(xlim_left, xlim_right)
    ax.set_ylim(ylim_bot, ylim_top)
    ax.set_xlabel('Timeline (days)')
    ax.set_ylabel('Survival')
    plt.grid()
    plt.show()
    return ax


def get_tree_height(tree):
    nodes = list(tree.nodes())
    nodes.sort(reverse=True)
    pos = 0  # 0th position is a recursion level
    return int(nodes[0][pos])


def make_layout(groups, node_size, root=dsc.root):
    def calc_pos_current_node(pos_parent, level, node_sz, sd):
        # sd == -1 if left; == 1 if right
        def calc_step_down(lvl, nd_sz):
            return pow(2., lvl)

        def calc_step_side(lvl, nd_sz):
            return pow(2., lvl)

        step_down = calc_step_down(level, node_sz)
        step_side = calc_step_side(level, node_sz)
        x_coord = pos_parent[0] - step_side if sd < 0 else pos_parent[0] + step_side
        return x_coord, pos_parent[1] - step_down

    def node_location(nd):
        return -1 if nd.split('-')[2] == 'left' else 1

    def node_level(nd):
        return int(nd.split('-')[0])

    pos = {root: (0, 0)}
    q = Queue()
    q.put(root)
    height = get_tree_height(groups)
    while not q.empty():
        node = q.get()
        scs = groups.successors(node)
        for nod in scs:
            pos[nod] = calc_pos_current_node(pos[node], height + 1 - node_level(nod),
                                             node_sz=node_size, sd=node_location(nod))
            q.put(nod)
    return pos


def make_edge_labels(groups):
    e = groups.edges()
    labels = [((u, v), get_edge_label(groups, u, v)) for u, v in e]
    # labels = [(e[i], values[i]) for i in xrange(len(e))]
    return {key: val for (key, val) in labels}


def plot_bin_tree(groups, node_size=50, with_labels=False, node_shape='o'):
    plt.figure()
    ax = plt.subplot(111)
    pos = make_layout(groups, node_size=node_size, root=dsc.root)
    # edge_labels = make_edge_labelgit s(groups)
    draw_networkx_nodes(groups, pos=pos, ax=ax, with_labels=with_labels, node_size=node_size, node_shape=node_shape,
                        node_color='k')
    draw_networkx_edges(groups, pos=pos, ax=ax, edge_color='k', arrows=False)
    # draw_networkx_edge_labels(groups, pos=pos, ax=ax, edge_labels=edge_labels, edge_color='k', rotate=True)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()
    return ax
