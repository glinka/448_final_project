import calcVotingModelProps as cvmp
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import matplotlib.pyplot as plt
import os

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    f = open(path_to_file, "r")
    params_str = f.readline()
    params = get_header_data(params_str)
    f.close()
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    print params
    return data, params

def get_header_data(header_str):
    BEGIN = 0
    comma = 1
    #create dict from header, based on key=value format in csv
    params = {}
    while comma > 0:
        equals = header_str.find("=")
        comma = header_str.find(",")
        params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
        header_str = header_str[comma+1:]
    params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
    #make integer, may not work especially well
    for key in params:
        if(params[key] % 1 == 0):
            params[key] = int(params[key])
    return params

def make_folder(base_name):
    folder_count = 0
    new_folder = base_name + str(folder_count) + '/'
    while os.path.exists(new_folder):
        folder_count = folder_count + 1
        new_folder = base_name + str(folder_count) + '/'
    os.mkdir(new_folder)
    return new_folder

def make_filename(base_name, params, unique_id=''):
    filename = base_name
    for key in params.keys():
        filename = filename + '_' + key + '_' + str(params[key])
    if not unique_id:
        filename = filename + '_' + unique_id
    return filename
    
def plot_timecourse(adj_data, opns_data, time_data, nvms_data, params, scalar_fn, average=False, ax=''): 
    n = params['n']
    n_data = time_data.shape[0]
    n_vms = params['nVms']
    starting_nvms = params['nVms']
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=n_vms)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Simulation step')
        ax.set_ylabel('Conflicts')
    ydata = []
    for i in range(n_data):
        #temp fp fix
        #n_vms = int(nvms_data[i] + 0.5)
        adj_matrices = [adj_data[(n_vms*i+k)*n:(n_vms*i+k+1)*n,:] for k in range(n_vms)]
        opn_vectors = [opns_data[(n_vms*i+k)*n:(n_vms*i+k+1)*n] for k in range(n_vms)]
        fn_evals = [scalar_fn(adj_matrices[k], opn_vectors[k]) for k in range(n_vms)]
        if average:
            ydata.append(np.average(fn_evals))
        else:
            ydata.append(fn_evals)
    if average:
        ax.plot(time_data, ydata)
    else:
        ydata = np.array(ydata)
        for k in range(starting_nvms):
            ax.plot(time_data, ydata[:,k], color=scalarMap.to_rgba(k))
    plt.show()

def plot_conflicts(data, params):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Minority fraction')
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.set_ylabel('Conflicts')
    ax2.set_xlabel('Simulation step')
    n = 1.0*params['n']
    ax1.plot(data[:,0], data[:,1]/n)
    ax2.plot(data[:,0], data[:,2])
    plt.show()

def plot_phase(data, params):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Minority fraction')
    ax.set_ylabel('Conflicts')
    n = 1.0*params['n']
    ax.plot(data[:,1]/n, data[:,2])
    plt.show()

def calc_conflict_triangles(adj, opns):
    n = adj.shape[0]
    tip_opn = 1
    ntris = 0
    for i in range(n):
        if(opns[i] == tip_opn):
            for j in range(n):
                if((adj[i][j] != 0) & (opns[j] != opns[i])):
                    for k in range(n):
                        if((adj[k][j] != 0) & (opns[k] == opns[j])):
                            if(adj[k][i] != 0):
                                ntris = ntris + 1
    return ntris/2

def plot_tc(adjs, opns, time, params):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    nentries = len(opns)
    n = params['n']
    ndata = nentries/n
    to_plot = []
    for i in range(ndata):
        adj = adjs[i*n:(i+1)*n,:]
        to_plot.append(calc_conflict_triangles(adj, opns[i*n:(i+1)*n]))
    ax.plot(time, to_plot)
    plt.show()

def plot_bifdiag(data, params):
    markers = ["p", "v", "D", "x", "s", "+", "*", "8", ".", "^"]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Probability of reattachment, ' + r"$\alpha$")
    ax.set_xlim((0,1))
    ax.set_ylabel('Final minority fraction')
    ax.set_ylim((0,0.5))
    ax.hold(True)
    n = data.shape[0]
    nPerSet = 11
    nData = int(n/nPerSet)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=nPerSet)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    for j in range(nData):
        i = nData-j-1
        toPlot = data[i*nPerSet:(i+1)*nPerSet,:]
        toPlot = np.array(filter(lambda x: (x[2] == 0), toPlot))
        print toPlot
        ax.scatter(toPlot[:,0], toPlot[:,1], color='k', marker=markers[i], label='Initial fraction: ' + str(data[i*nPerSet,3]))
    ax.legend(loc=2)
    plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-cpic', '--plot-cpi-conflicts', action='store_true', default=False)
    parser.add_argument('-cpimf', '--plot-cpi-minority-fraction', action='store_true', default=False)
    parser.add_argument('-cpipp', '--plot-cpi-phase-portrait', action='store_true', default=False)
    parser.add_argument('-avg', '--average', action='store_true', default=False)
    parser.add_argument('-tc', '--plot-tc', action='store_true', default=False)
    parser.add_argument('-phase', '--plot-phase-space', action='store_true', default=False)
    parser.add_argument('-bif', '--plot-bif-diag', action='store_true', default=False)
    parser.add_argument('-tris', '--plot-triangles', action='store_true', default=False)
    args = parser.parse_args()
    #should have one of each file: CPIAdj_... and CPIOpns_...
    #which will have the same header (and thus same params)
    for filename in args.input_files:
        if 'Adj' in filename:
            adj_data, params = get_data(filename)
        elif 'Opns' in filename:
            opns_data, params = get_data(filename)
        elif 'Times' in filename:
            time_data, params = get_data(filename)
        elif 'ids' in filename:
            nvms_data, params = get_data(filename)
        elif 'graphStats' in filename:
            graphStats, params = get_data(filename)
        elif 'bifData' in filename:
            bifdata, params = get_data(filename, header_rows=0)
    if args.plot_cpi_conflicts:
        plot_timecourse(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_conflicts, average=args.average)
    if args.plot_cpi_minority_fraction:
        plot_timecourse(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_minority_fraction, average=args.average)
    if args.plot_cpi_phase_portrait:
        plot_phase_portrait(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_minority_fraction, cvmp.get_conflicts)
    if args.plot_tc:
        plot_conflicts(graphStats, params)
    if args.plot_phase_space:
        plot_phase(graphStats, params)
    if args.plot_bif_diag:
        plot_bifdiag(bifdata, params)
    if args.plot_triangles:
        plot_tc(adj_data, opns_data, time_data, params)
    testadj = np.array([[0,1,1],[1,0,1],[1,1,0]])
    testopns = np.array([1,2,2])
    print calc_conflict_triangles(testadj, testopns)
