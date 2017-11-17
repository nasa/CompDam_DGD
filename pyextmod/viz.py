import numpy as np
from mpl_toolkits.mplot3d import Axes3D, proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
import re

class Points():

    def __init__(self, lc, alpha, F, F_bulk):

        self.F = F
        self.F_bulk = F_bulk

        self.initial = {}
        self.deformed = {}

        lc1 = lc[0]
        lc2 = lc[1]
        lc3 = lc[2]

        # Initial configuration of the total volume
        self.initial['total'] = np.array([
                [0,   0,   0],
                [lc1, 0,   0],
                [lc1, lc2, 0],
                [0,   lc2, 0],
                [0,   0,   lc3],
                [lc1, 0,   lc3],
                [lc1, lc2, lc3],
                [0,   lc2, lc3]
            ])

        # Initial configuration of one half of the bulk
        self.initial['b1'] = np.array([
                [0,   0,                         0],
                [lc1, 0,                         0],
                [lc1, lc2/2+lc3/2*np.tan(alpha), 0],
                [0,   lc2/2+lc3/2*np.tan(alpha), 0],
                [0,   0,                         lc3],
                [lc1, 0,                         lc3],
                [lc1, lc2/2-lc3/2*np.tan(alpha), lc3],
                [0,   lc2/2-lc3/2*np.tan(alpha), lc3]
            ])

        # Initial configuration of the other half of the bulk
        self.initial['b2'] = np.array([
                [0,   lc2/2+lc3/2*np.tan(alpha), 0],
                [lc1, lc2/2+lc3/2*np.tan(alpha), 0],
                [lc1, lc2, 0],
                [0,   lc2, 0],
                [0,   lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2, lc3],
                [0,   lc2, lc3]
            ])

        # Initial configuration of the cohesive element
        self.initial['coh'] = np.array([
                [0,   lc2/2+lc3/2*np.tan(alpha), 0],
                [0,   lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2/2+lc3/2*np.tan(alpha), 0],
                [0,   lc2/2+lc3/2*np.tan(alpha), 0],
                [0,   lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2/2-lc3/2*np.tan(alpha), lc3],
                [lc1, lc2/2+lc3/2*np.tan(alpha), 0]
            ])

        self.deformed['total'] = np.transpose(np.matmul(F, np.transpose(self.initial['total'])))
        self.deformed['b1'] = np.transpose(np.matmul(F_bulk, np.transpose(self.initial['b1'])))
        self.deformed['b2'] = np.transpose(np.matmul(F_bulk, np.transpose(self.initial['b2'])))
        self.deformed['b2'][:,1] = self.deformed['b2'][:,1] + self.deformed['total'][3,1] - self.deformed['b2'][3,1]
        self.deformed['b2'][:,2] = self.deformed['b2'][:,2] + self.deformed['total'][3,2] - self.deformed['b2'][3,2]

        self.deformed['coh'] = np.array([
                self.deformed['b1'][3,:],
                self.deformed['b1'][7,:],
                self.deformed['b1'][6,:],
                self.deformed['b1'][2,:],
                self.deformed['b2'][0,:],
                self.deformed['b2'][4,:],
                self.deformed['b2'][5,:],
                self.deformed['b2'][1,:]
        ])


def visualize(lc, alpha, F, pathToLogFile, initMD=1, initEQk=1):
    # Load the data
    d = _load_data(lc, alpha, F, pathToLogFile)
    md_vals = d.keys()

    # Initialize the figure
    fig = plt.figure(figsize=(15,10))
    ax1 = plt.subplot2grid((12, 12), (0, 0), rowspan=8, colspan=6, projection='3d')
    ax2 = plt.subplot2grid((12, 12), (0, 7), rowspan=7, colspan=6)
    ax_md_label = plt.subplot2grid((12, 12), (9, 1))
    ax_md_goto = plt.subplot2grid((12, 12), (9, 3))
    ax_md_prev = plt.subplot2grid((12, 12), (9, 4))
    ax_md_next = plt.subplot2grid((12, 12), (9, 5))
    ax_eq_label = plt.subplot2grid((12, 12), (10, 1))
    ax_eq_goto = plt.subplot2grid((12, 12), (10, 3))
    ax_eq_prev = plt.subplot2grid((12, 12), (10, 4))
    ax_eq_next = plt.subplot2grid((12, 12), (10, 5))
    ax_view12 = plt.subplot2grid((12, 12), (11, 1))
    ax_view23 = plt.subplot2grid((12, 12), (11, 2))
    ax_view13 = plt.subplot2grid((12, 12), (11, 3))
    plt.subplots_adjust(wspace=0.4)

    # Generate the initial visualization
    ax1 = visualize_iteration(ax1, d[initMD][initEQk]['pts'])
    ax1.grid(False)

    # Plot error vs EQk
    def plot_err(md):
        err_data = [v['err'] for k,v in d[md].iteritems()]
        ax2.semilogy(err_data, '-k+')
        ax2.set_xlabel('EQk')
        ax2.set_ylabel('Err')
    plot_err(initMD)

    # Calculate max ranges
    mdmax = max(d.keys())
    eqkmax = max(d[initMD].keys())

    # Text boxes
    ax_md_label.axis('off')
    txt_md = ax_md_label.text(0.5, 0.5, 'MD:  ' + str(initMD) + ' of ' + str(mdmax), fontsize=12)
    ax_eq_label.axis('off')
    txt_eq = ax_eq_label.text(0.5, 0.5, 'EQk: ' + str(initEQk) + ' of ' + str(eqkmax), fontsize=12)

    # Frame change buttons event handlers
    class frame():

        def __init__(self, mdinit, mdmax, eqkinit, eqkmax):
            self.md = mdinit
            self.mdmax = mdmax
            self.eqk = eqkinit
            self.eqkmax = eqkmax

        def next_md(self, event):
            if self.md < self.mdmax:
                self.md += 1
                self.eqk = 1
                self.eqkmax = max(d[self.md].keys())
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                ax2.clear()
                plot_err(self.md)
                txt_md.set_text('MD:  ' + str(self.md) + ' of ' + str(self.mdmax))
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()

        def prev_md(self, event):
            if self.md > 1:
                self.md -= 1
                self.eqk = 1
                self.eqkmax = max(d[self.md].keys())
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                ax2.clear()
                plot_err(self.md)
                txt_md.set_text('MD:  ' + str(self.md) + ' of ' + str(self.mdmax))
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()

        def goto_md(self, text):
            i = int(text)
            if i >= 1 and i <= self.mdmax:
                self.md = i
                self.eqk = 1
                self.eqkmax = max(d[self.md].keys())
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                ax2.clear()
                plot_err(self.md)
                txt_md.set_text('MD:  ' + str(self.md) + ' of ' + str(self.mdmax))
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()

        def next_eqk(self, event):
            if self.eqk < self.eqkmax:
                self.eqk += 1
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()

        def prev_eqk(self, event):
            if self.eqk > 1:
                self.eqk -= 1
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()

        def goto_eqk(self, text):
            i = int(text)
            if i >= 1 and i <= self.eqkmax:
                self.eqk = i
                ax1.clear()
                visualize_iteration(ax1, d[self.md][self.eqk]['pts'])
                txt_eq.set_text('EQk: ' + str(self.eqk) + ' of ' + str(self.eqkmax))
                fig.canvas.draw_idle()


    # Buttons for changing MD and EQk
    callback = frame(mdinit=initMD, mdmax=mdmax, eqkinit=initEQk, eqkmax=eqkmax)
    bnext_md = Button(ax_md_next, 'Next')
    bnext_md.on_clicked(callback.next_md)
    bprev_md = Button(ax_md_prev, 'Previous')
    bprev_md.on_clicked(callback.prev_md)
    bnext_eqk = Button(ax_eq_next, 'Next')
    bnext_eqk.on_clicked(callback.next_eqk)
    bprev_eqk = Button(ax_eq_prev, 'Previous')
    bprev_eqk.on_clicked(callback.prev_eqk)

    tb_md = TextBox(ax_md_goto, 'Go to ', initial='MD')
    tb_md.on_submit(callback.goto_md)
    tb_eq = TextBox(ax_eq_goto, 'Go to ', initial='EQk')
    tb_eq.on_submit(callback.goto_eqk)

    # Buttons for changing view
    def view12(event):
        ax1.view_init(elev=90., azim=-90.)
    def view23(event):
        ax1.view_init(elev=0., azim=0.)
    def view13(event):
        ax1.view_init(elev=0., azim=-90.)
    bview12 = Button(ax_view12, 'View 1-2')
    bview12.on_clicked(view12)
    bview23 = Button(ax_view23, 'View 2-3')
    bview23.on_clicked(view23)
    bview13 = Button(ax_view13, 'View 1-3')
    bview13.on_clicked(view13)

    plt.show()


def visualize_iteration(ax, points):

    _plot_parallelpiped(ax, points.deformed['total'], alpha=0, edgecolor='k')
    # plot_parallelpiped(ax, deformed_pts, edgecolor='c', alpha=0)
    _plot_parallelpiped(ax, points.deformed['b1'], facecolor=[0.5, 0.5, 1], edgecolor='r')
    _plot_parallelpiped(ax, points.deformed['b2'], facecolor=[0.5, 0.5, 1], edgecolor='b')
    _plot_parallelpiped(ax, points.deformed['coh'], facecolor='xkcd:yellow', edgecolor='k', ls='dotted')

    # Axis label
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    return ax


def _load_data(lc, alpha, F, pathToLogFile):

    # Store the output as a nest list
    data = {}
    # regexs
    md = re.compile(r' MD = (\d+)')
    eqk = re.compile(r' EQk = (\d+)')
    err = re.compile(r' err = ([\d.E+-]+)')
    fb = re.compile(r' Fb =')

    # Load the log file
    with open(pathToLogFile, 'r') as lf:
        lines = lf.readlines()

    # Parse the log file
    numLines = len(lines)
    i = 0
    while (i < numLines):
        line = lines[i]

        # MD
        md_match = md.match(line)
        if md_match:
            current_md = int(md_match.groups()[0])
            i += 1
            continue

        # EQk
        eqk_match = eqk.match(line)
        if eqk_match:
            current_eqk = int(eqk_match.groups()[0])
            i += 1
            continue

        # err
        err_match = err.match(line)
        if err_match:
            current_err = float(err_match.groups()[0])
            i += 1
            continue

        # Fb
        fb_match = fb.match(line)
        if fb_match:
            F_bulk = np.array([
                [float(x) for x in lines[i+1].split(',')],
                [float(x) for x in lines[i+2].split(',')],
                [float(x) for x in lines[i+3].split(',')]
            ])
            if current_md not in data:
                data[current_md] = {}
            data[current_md][current_eqk] = {'pts': Points(lc, alpha, F, F_bulk), 'err': current_err}
            i += 4
            continue

        # No matches
        i += 1

    return data


def _plot_parallelpiped(ax, points, facecolor=[0.5, 0.5, 1], edgecolor='k', alpha=0.2, ls='solid'):
    """
    Adds transparent parallelpiped in 3d to the axes ax
    points is n x 3 array of vertex coordinates
    """

    # list of sides' polygons of figure
    verts = [[points[0],points[1],points[2],points[3]],
        [points[4],points[5],points[6],points[7]],
        [points[0],points[1],points[5],points[4]],
        [points[2],points[3],points[7],points[6]],
        [points[1],points[2],points[6],points[5]],
        [points[4],points[7],points[3],points[0]],
        [points[2],points[3],points[7],points[6]]]

    # plot vertices
    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2])

    # plot sides
    collection = Poly3DCollection(verts, linewidths=1, edgecolors=edgecolor, alpha=alpha)
    collection.set_linestyle(ls)
    collection.set_facecolor(facecolor)
    ax.add_collection3d(collection)
