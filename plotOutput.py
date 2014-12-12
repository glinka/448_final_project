import calcVotingModelProps as cvmp
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    f = open(path_to_file, "r")
    params_str = f.readline()
    params = get_header_data(params_str)
    f.close()
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    return data, params

def read_ids(filename, header_rows=1, conversion=int, **kwargs):
    path_to_file = os.path.realpath(filename)
    f = open(path_to_file, "r")
    params_str = ''
    params = {}
    str_ids = []
    [str_ids.append(line) for line in f]
    if header_rows == 1:
        params_str = str_ids[0]
        params = get_header_data(params_str)
        str_ids = str_ids[header_rows:]
    f.close()
    ids = []
    string_numbers = [str(i) for i in range(10)]
    for line in str_ids:
        myline = line
        ids.append([])
        comma_loc = 0
        while comma_loc >= 0:
            comma_loc = myline.find(",")
            ids[-1].append(conversion(myline[0:comma_loc]))
            myline = myline[comma_loc+1:]
    return np.array(ids), params

def thin_array(array, frac_to_keep=0.5, new_npts=None):
    # will keep 100*frac_to_keep% of the data 
    # in each column of array, evenly spaced
    # thus, array has form col1 col2 col3 .. colN
    # and returns col1 col2...colN reduced to frac_to_keep of
    # original length
    # no thinning is possible for frac_to_keep > 0.5,
    # at least in this arrangement
    if new_npts is None:
        if frac_to_keep > 0.5:
            return array

        npts = array.shape[0]
        new_npts = int(frac_to_keep*npts)
        spacing = npts/new_npts
        ncols = array.shape[1]
        thinned_array = np.zeros((new_npts, ncols))
    else:
        npts = array.shape[0]
        if new_npts > npts/2:
            return array
        npts = array.shape[0]
        spacing = npts/new_npts
        ncols = array.shape[1]
        thinned_array = np.zeros((new_npts, ncols))
    for j in range(ncols):
        for i in range(new_npts):
            thinned_array[i,j] = array[spacing*i, j]
    return thinned_array

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

def plot_timecourse(adj_data, opns_data, time_data, ids_data, params, folder, scalar_fn, average=False, ax=''):
    starting_nvms = params['nVms']
    n = params['n']
    n_data = time_data.shape[0]
    ydata = {}
    for i in range(starting_nvms):
        ydata[str(i)] = []
    data_count = 0
    for i in range(n_data):
        for j in ids_data[i]:
            adj = adj_data[data_count*n:(data_count+1)*n,:]
            opns = opns_data[data_count*n:(data_count+1)*n]
            ydata[str(j)].append(scalar_fn(adj, opns))
            data_count = data_count + 1
    final_ydata = []
    if average:
        for i in range(n_data):
            final_ydata.append(np.average([ydata[str(j)][i] for j in ids_data[i]]))
    else:
        for key in ydata.keys():
            final_ydata.append(ydata[key])
    final_ydata = np.array(final_ydata)
    if not ax:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    if average:
        ax.plot(time_data, final_ydata)
    else:
        ax.hold(True)
        for i in range(starting_nvms):
            npoints = len(final_ydata[i])
            ax.plot(time_data[:npoints], final_ydata[i])
#    plt.show()

def plot_conflicts(conflicts_data, time_data, ids_data, params, folder='', average=False, ax=''):
    starting_nvms = params['nVms']
    n = params['n']
    n_data = time_data.shape[0]
    ydata = {}
    for i in range(starting_nvms):
        ydata[str(i)] = []
    data_count = 0
    for i in range(n_data):
        count = 0
        for j in ids_data[i]:
            ydata[str(j)].append(conflicts_data[i][count])
            count = count + 1
    final_ydata = []
    if average:
        for i in range(n_data):
            to_average = []
            [to_average.append(ydata[str(j)][i]) for j in ids_data[i]]
            #pad with zeros for finished runs
            [to_average.append(0) for j in range(starting_nvms - len(ids_data[i]))]
            final_ydata.append(np.average(to_average))
    else:
        for key in ydata.keys():
            final_ydata.append(ydata[key])
    final_ydata = np.array(final_ydata)
    if not ax:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    if average:
        # colors jumps in red
        ntimes = time_data.shape[0]
        to_plot = []
        jump_size = 10000
        lastindex = 0
        for i in range(ntimes-1):
            if (time_data[i+1] - time_data[i]) > jump_size:
                ax.plot(time_data[lastindex:i], final_ydata[lastindex:i], lw=2, c='b')
                ax.plot(time_data[i:i+2], final_ydata[i:i+2], lw=3, c='r')
                lastindex=i+2
        ax.plot(time_data[lastindex:], final_ydata[lastindex:], label='CPI', c='b', lw=2)
        # normally use this:
        # ax.plot(time_data, final_ydata, label='Cpi', lw=2)
    else:
        ax.hold(True)
        for i in range(starting_nvms):
            npoints = len(final_ydata[i])
            ax.plot(time_data[:npoints], final_ydata[i])
    ax.legend()
#    plt.show()

def plot_minorities(minorities_data, time_data, params, folder='', average=True, ax=''):
    #we'll assume we can simply import the minorities data w/ genfromtext
    starting_nvms = params['nVms']
    n = params['n']
    n_data = time_data.shape[0]
    ydata = {}
    final_ydata = []
    if average:
        final_ydata = np.average(minorities_data, 1)
    else:
        final_ydata = np.array(minorities_data)
    if not ax:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    if average:
        # specifically designed to color jumps differently
        ntimes = time_data.shape[0]
        to_plot = []
        jump_size = 10000
        lastindex = 0
        for i in range(ntimes-1):
            if (time_data[i+1] - time_data[i]) > jump_size:
                ax.plot(time_data[lastindex:i], final_ydata[lastindex:i], lw=2, c='b')
                ax.plot(time_data[i:i+2], final_ydata[i:i+2], lw=3, c='r')
                lastindex=i+2
        ax.plot(time_data[lastindex:], final_ydata[lastindex:], label='CPI', c='b', lw=2)
        # use this normally:  
        # ax.plot(time_data, final_ydata, label='Cpi', lw=2)
    else:
        ax.hold(True)
        for i in range(starting_nvms):
            ax.plot(time_data[:], final_ydata[:,i])
##############################
##############################
    #run specific, only for use on csv_data_fast12_/*
##############################
##############################
    time_start_enlarge = 35000
    time_end_enlarge = 60000
    data_to_enlarge = np.zeros([n_data, 2])
    data_to_enlarge[:,0] = final_ydata
    data_to_enlarge[:,1] = time_data
    data_to_enlarge.shape = [n_data, 2]
    #filter out times outside of desired interval
    data_to_enlarge = filter(lambda t: ((t[1] >= time_start_enlarge) & (t[1] <= time_end_enlarge)), data_to_enlarge)
    data_to_enlarge = np.array(data_to_enlarge)
    n_enlarged_data = data_to_enlarge.shape[0]
    minorities_to_enlarge = data_to_enlarge[:,0]
    minor_frac_ax = (np.max(minorities_to_enlarge) - np.min(minorities_to_enlarge))/0.5
    times_to_enlarge = data_to_enlarge[:,1]
    time_frac_ax = (np.max(times_to_enlarge) - np.min(times_to_enlarge))/250000.
    scale_factor = 0.2/time_frac_ax
    fig = plt.gcf()
    embedded_ax = fig.add_axes([0.5, 0.4, 0.18, scale_factor*minor_frac_ax])
    A = np.ones([n_enlarged_data, 2])
    A[:,0] = times_to_enlarge
    rhs = minorities_to_enlarge
    m, b = np.linalg.lstsq(A, rhs)[0]
    print m, b
    embedded_ax.set_xticks([])
    embedded_ax.set_yticks([])
    embedded_ax.plot(times_to_enlarge, minorities_to_enlarge)
    embedded_ax.plot(times_to_enlarge, m*times_to_enlarge+b, c='r', label='linear fit')
    # embedded_ax.legend()
    example_times = np.linspace(0, 250000, 100)
    # ax.plot(example_times, -2.7e-6*example_times+0.379558, label='1')
    # ax.plot(example_times, -2.81e-7*example_times+0.236155, label='2')
    # ax.plot(example_times, -2.342e-6*example_times+0.45, label='3')
    # ax.plot(example_times, -5.85e-7*example_times+0.23, label='4')
    # ax.set_ylim([0,0.5])
    ax.legend()
    # plt.show()

def plot_single_conflicts(data, params, ax='', average=True):
    # textsize = 22
    # ax1 = fig.add_subplot(211)
    # ax1.set_ylabel('Minority fraction', fontsize=textsize)
    # plt.tick_params(axis='both', which='major', labelsize=18)
    # ax1.set_ylim((0,0.5))
    # ax1.set_yticks([(x+1)/10.0 for x in range(5)])
    # ax2 = fig.add_subplot(212, sharex=ax1)
    # ax2.set_ylabel('Conflicts', fontsize=textsize)
    # ax2.set_ylim((0,400))
    # ax2.set_yticks([(x+1)*100 for x in range(4)])
    # ax2.set_xlabel('Simulation step', fontsize=textsize)
    # plt.tick_params(axis='both', which='major', labelsize=18)
    # n = 1.0*params['n']
    # ax1.plot(data[:,0], data[:,1]/n)
    # ax2.plot(data[:,0], data[:,2])
    # plt.show()
    toplot = []
    nruns = len(data)
    maxiter = np.max([data[j][0].shape[0] for j in range(nruns)])
    maxstep = int(np.max([data[j][0][-1,0] for j in range(nruns)]))
    runlengths = []
    for i in range(nruns):
        runlengths.append(data[i][0].shape[0])
    if not ax:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        textsize = 22
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.xticks(rotation=20)
        ax.set_xlabel('simulation step', fontsize=textsize)
        ax.set_ylabel('conflicts', fontsize=textsize)
    n = 1.0*params['n']
    if average:
        for i in range(maxiter):
            to_average = []
            for j in range(nruns):
                if runlengths[j] > i:
                    to_average.append(data[j][0][i,2])
                else:
                    to_average.append(0)
            toplot.append(np.average(to_average))
        ax.plot(np.linspace(0, maxstep, maxiter), toplot, label='Averaged direct simulations', lw=2, c='g')
    else:
        for j in range(nruns):
            ax.plot(data[j][0][:,0], data[j][0][:,2])



def plot_single_minorities(data, params, ax='', average=True):
    # textsize = 22
    # ax1 = fig.add_subplot(211)
    # ax1.set_ylabel('Minority fraction', fontsize=textsize)
    # plt.tick_params(axis='both', which='major', labelsize=18)
    # ax1.set_ylim((0,0.5))
    # ax1.set_yticks([(x+1)/10.0 for x in range(5)])
    # ax2 = fig.add_subplot(212, sharex=ax1)
    # ax2.set_ylabel('Conflicts', fontsize=textsize)
    # ax2.set_ylim((0,400))
    # ax2.set_yticks([(x+1)*100 for x in range(4)])
    # ax2.set_xlabel('Simulation step', fontsize=textsize)
    # plt.tick_params(axis='both', which='major', labelsize=18)
    # n = 1.0*params['n']
    # ax1.plot(data[:,0], data[:,1]/n)
    # ax2.plot(data[:,0], data[:,2])
    # plt.show()
    toplot = []
    nruns = len(data)
    maxiter = np.max([data[j][0].shape[0] for j in range(nruns)])
    maxstep = int(np.max([data[j][0][-1,0] for j in range(nruns)]))
    runlengths = []
    for i in range(nruns):
        runlengths.append(data[i][0].shape[0])
    if not ax:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    ax.set_ylabel('Minorities', fontsize=textsize)
    ax.set_xlabel('Iterations', fontsize=textsize)
    n = 1.0*params['n']
    if average:
        for i in range(maxiter):
            to_average = []
            for j in range(nruns):
                if runlengths[j] > i:
                    to_average.append(data[j][0][i,1]/n)
                else:
                    to_average.append(data[j][0][-1,1]/n)
            toplot.append(np.average(to_average))
        ax.plot(np.linspace(0, maxstep, maxiter), toplot, label='Averaged direct simulations', lw=2, c='g')
    else:
        for j in range(nruns):
            ax.plot(data[j][0][:,0], data[j][0][:,1])

def plot_phase_portrait(data, params):
    ndata = len(data)
    n = params['n']
    # make into fractions instead of counts
    for run in data:
        run[:,1] = (1.0*run[:,1])/n
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    textsize = 36
    labelsize = 30
    ax.set_xlabel('Fraction voters with first opinion', fontsize=textsize)
    ax.set_ylabel('Conflicts', fontsize=textsize)
    ax.set_color_cycle([[1-i, np.cos(np.pi*i/2), i] for i in np.array(range(ndata))/float(ndata)])
    for run in data:
        ax.plot(run[:,1], run[:,2])
    plt.tick_params(axis='both', which='major', labelsize=labelsize)
    plt.show()

def plot_cpi_phase_portrait(minorities, conflicts, params, avg=True, ax=None):
    if ax is None:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    ax.set_xlabel('Minorities')
    ax.set_ylabel('Conflicts')
    ax.plot([np.average(minorities_slice) for minorities_slice in minorities], [np.average(conflict_slice) for conflict_slice in conflicts])
    plt.show()

def plot_triangles(data, params):
    n = 1.0*params['n']
    triangles = data[:,3]
    minorities = data[:,1]/n
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.scatter(minorities, triangles, lw=0, alpha=0.7)
    plt.show(fig)

def plot_simple_conflicts(data, params):
    n = 1.0*params['n']
    triangles = data[:,3]
    conflicts = data[:,2]
    minorities = data[:,1]/n
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.scatter(minorities, conflicts, lw=0, alpha=0.7)
    plt.show(fig)

def progress_bar(current_iter, total_iter, elapsed_time=None):
    perc = int((100.0*current_iter)/total_iter)
    percf = (100.0*current_iter)/total_iter
    bar = '\r['
    for i in range(perc):
        bar = bar + '|'
    for i in range(100-perc):
        bar = bar + ' '
    bar = bar + '] '
    if elapsed_time is not None:
        bar = bar + str(int(elapsed_time/(percf/100.0)) - int(elapsed_time)) + 's remaining'
    print bar,
    sys.stdout.flush()


def plot_for_anim(data, params, npts=None, thinning=1, filemarker=''):
    # reduce number of pts used if desired
    if npts is not None:
        thinned_data = []
        for run in data:
            thinned_data.append(thin_array(run, new_npts=npts))
        data = thinned_data
    else:
        thinned_data = []
        for run in data:
            thinned_data.append(thin_array(run, frac_to_keep=thinning))
        data = thinned_data
    nruns = len(data)
    n = 1.0*params['n']
    nsteps = [data[i].shape[0] for i in range(nruns)]
    for run in data:
        run[:,1] = run[:,1]/n
    # set up figure for plotting
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    textsize = 20
    lablesize = 4
    ax.set_xlabel('Fraction voters with first opinion', fontsize=textsize)
    ax.set_ylabel('Conflicts', fontsize=textsize)
    ax.set_xlim((0,1))
    ax.set_ylim((0, 1.1*[np.max(run[:,2]) for run in data][0]))
    # plt.tick_params(axis='both', which='major', labelsize=labelsize)
    fadeFromWhiteCmaps = ['Blues', 'YlOrRd', 'BuGn', 'BuPu', 'YlGn', 'Greens', 'RdPu', 'OrRd', 'Greys', 'Purples']
    ncmaps = len(fadeFromWhiteCmaps)
    start = time.clock()
    maxsteps = np.max(nsteps)
    endtimes = [run[-1,0] for run in data]
    # plot every step in a new image,
    # if one of the runs has completed
    # fade it to white
    for i in range(1,maxsteps):
        j = 0
        for run in data:
            s = []
            if i < nsteps[j]:
                ax.scatter(run[:i,1], run[:i,2], c=run[:i,0], cmap=fadeFromWhiteCmaps[j], lw=0, alpha=0.7, label='traj ' + str(j))
            else:
                fade = (1.0*i)/nsteps[j]
                ax.scatter(run[:,1], run[:,2], c=run[:,0]/fade, cmap=fadeFromWhiteCmaps[j], vmax=endtimes[j], lw=0, alpha=0.7, label='traj ' + str(j))
            j = (j + 1) % ncmaps
            ax.legend()
        plt.savefig('tempanim/' + filemarker + 'pp_' + str(i) + '.png')
        progress_bar(i+1, maxsteps, time.clock() - start)
    print
        # okay to remove elements while looping?
        # del s
        # for j in range(len(s)):
        #     s[j].remove()
        #     del s[j]


def vm_dmaps_embedding(eigvals, eigvects, params, t=0, plot_2d=True, plot_3d=False, label=''):
    from mpl_toolkits.mplot3d import Axes3D
    ntypes = 1 # params['ninits']
    n = eigvects.shape[1]
    n_pertype = n/ntypes
    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    eigvects_to_plot = np.array([eigvects[-i,:] for i in range(2, eigvects.shape[0] + 1)])
    eigvals_to_plot = np.array([eigvals[-i] for i in range(2, eigvals.shape[0] + 1)])
    nvects = eigvals.shape[0] - 1
    output_filename = "vm_embedding_" + label
    cmaps = ['jet', 'autumn', 'spectral', 'winter', 'summer', 'PrOr', 'RdBu', 'RdYlBu', 'RdYlGn']

    mfs = np.genfromtxt("./single_runs/mf_evos.csv")

    if plot_2d is False:
        print 'Not plotting 2d embeddings'
    else:
        print 'Plotting 2d embeddings'
        # plot 2d embeddings
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        n = eigvects.shape[1]
        for i in range(nvects):
            for j in range(i+1, nvects):
                xvals = np.power(eigvals_to_plot[i], t)*eigvects_to_plot[i,:]
                yvals = np.power(eigvals_to_plot[j], t)*eigvects_to_plot[j,:]
                for k in range(ntypes):
                    ax.scatter(xvals[k*n_pertype:(k+1)*n_pertype], yvals[k*n_pertype:(k+1)*n_pertype], c=mfs[k*n_pertype:(k+1)*n_pertype], lw=0, alpha=0.7, cmap=cmaps[k])# np.arange(n_pertype), lw=0, alpha=0.7, cmap=cmaps[k])
                    ax.hold(True)
                ax.set_xlim((np.min(xvals), np.max(xvals)))
                ax.set_ylim((np.min(yvals), np.max(yvals)))
                ax.set_xlabel('eigvect ' + str(i+1))
                ax.set_ylabel('eigvect ' + str(j+1))
                ax.set_title('vm dmaps embedding')
                plt.savefig("./figs/embeddings/dmaps/2d/" + output_filename + "eigvects_" + str(i+1) + str(j+1) + ".png")
                ax.hold(False)
    # plot 3d embeddings
    if plot_3d is False:
        print 'Not plotting 3d embeddings'
    else:
        print 'Plotting 3d embeddings'
        for i in range(nvects):
            for j in range(i+1, nvects):
                for k in range(j+1, nvects):
                    fig = plt.figure(facecolor='w')
                    ax = fig.add_subplot(111, projection='3d')
                    xvals = np.power(eigvals_to_plot[i], t)*eigvects_to_plot[i,:]
                    yvals = np.power(eigvals_to_plot[j], t)*eigvects_to_plot[j,:]
                    zvals = np.power(eigvals_to_plot[k], t)*eigvects_to_plot[k,:]
                    ax.hold(False)
                    # ax.scatter(xvals[:n/2], yvals[:n/2], zvals[:n/2], c=np.arange(n/2), lw=0, alpha=1, cmap='Reds')
                    # ax.hold(True)
                    # ax.scatter(xvals[-n/2:], yvals[-n/2:], zvals[-n/2:], c=np.arange(n/2), lw=0, alpha=1, cmap='jet')
                    ax.scatter(xvals, yvals, zvals, c=np.arange(n), lw=0, alpha=1, cmap='jet')
                    ax.set_xlim((np.min(xvals), np.max(xvals)))
                    ax.set_ylim((np.min(yvals), np.max(yvals)))
                    ax.set_zlim((np.min(zvals), np.max(zvals)))
                    ax.set_xlabel('eigvect ' + str(i+1))
                    ax.set_ylabel('eigvect ' + str(j+1))
                    ax.set_zlabel('eigvect ' + str(k+1))
                    plt.show(fig)
                    # plt.savefig("./figs/embeddings/dmaps/3d/" + output_filename + "eigvects_" + str(i+1) + str(j+1) + str(k+1) + ".png")

def conflicts_phaseplot(minority_fracs, conflicts, params):
    n = float(params['n'])
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.hold(True)
    ndatasets = len(minority_fracs)
    cs = ['b', 'g', 'r', 'c', 'b', 'g', 'r', 'c']
    for i in range(ndatasets):
        ax.plot(minority_fracs[i]/n, conflicts[i], c=cs[i])
        ax.set_xlabel('fraction opinion "a"')
        ax.set_ylabel('conflict squares')
    plt.show()

def plot_timecourses(minority_fracs, conflicts, cherry_conflicts, times, params):
    import matplotlib.gridspec as gs

    n = float(params['n'])
    fig = plt.figure(facecolor='w')
    gspec = gs.GridSpec(3,1)
    ax_top = fig.add_subplot(gspec[0])
    ax_mid = fig.add_subplot(gspec[1])
    ax_bottom = fig.add_subplot(gspec[2])
    ax_top.plot(times, cherry_conflicts)
    ax_mid.plot(times, conflicts)
    ax_bottom.plot(times, minority_fracs/n)
    plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-cpic', '--plot-cpi-conflicts', action='store_true', default=False)
    parser.add_argument('-sc', '--plot-single-conflicts', action='store_true', default=False)
    parser.add_argument('-cpim', '--cpi-minorities', action='store_true', default=False)
    parser.add_argument('-sm', '--single-minorities', action='store_true', default=False)
    parser.add_argument('-cpimf', '--plot-cpi-minority-fraction', action='store_true', default=False)
    parser.add_argument('-cpipp', '--plot-cpi-phase-portrait', action='store_true', default=False)
    parser.add_argument('-ppp', '--plot-phase-portrait', action='store_true', default=False)
    parser.add_argument('-pt', '--plot-triangles', action='store_true', default=False)
    parser.add_argument('--plot-animation', action='store_true', default=False)
    parser.add_argument('-avg', '--average', action='store_true', default=False)
    parser.add_argument('--dmaps-embeddings', action='store_true', default=False)
    parser.add_argument('--cherries', action='store_true', default=False)
    parser.add_argument('--squares', action='store_true', default=False)
    parser.add_argument('--timecourses', action='store_true', default=False)
    args = parser.parse_args()
    #should have one of each file: CPIAdj_... and CPIOpns_...
    #which will have the same header (and thus same params)
    graphdata = []
    # fig = plt.figure(facecolor='w')
    # myax = fig.add_subplot(111)
    # myax.hold(True)
    # textsize = 34
    # labelsize=30
    # myax.set_xlabel('Iterations', fontsize=textsize)
    # myax.set_ylabel('Conflicts', fontsize=textsize)
    # plt.xticks(rotation=20)
    # plt.tick_params(axis='both', which='major', labelsize=labelsize)
    # for filename in args.input_files:
    #     slash = filename.find('/')
    #     folder = filename[:slash+1]
    #     if 'Adj' in filename:
    #         adj_data, params = get_data(filename)
    #     if 'Opns' in filename:
    #         opns_data, cpiparams = get_data(filename)
    #     elif 'Times' in filename:
    #         time_data, params = get_data(filename)
    #     elif 'ids' in filename:
    #         ids_data, cpiparams = read_ids(filename)
    #     elif 'Conflicts' in filename:
    #         conflicts_data, cpiparams = read_ids(filename, header_rows=0)
    #     elif "graphstats" in filename:
    #         graphdata.append(get_data(filename))
    # if args.plot_cpi_conflicts:
    #     plot_timecourse(adj_data, opns_data, time_data, ids_data, params, folder, cvmp.get_conflicts, average=args.average)
    if args.plot_phase_portrait:
        data = []
        for f in args.input_files:
            d, params = get_data(f)
            data.append(d)
        plot_phase_portrait(data, params)
    if args.plot_animation:
        data = []
        for f in args.input_files:
            d, params = get_data(f)
            data.append(d)
        plot_for_anim(data, params, thinning=1.0/1000)
        # for special case of importing data with ci = 1
        short_data = []
        for run in data:
            short_data.append(run[:200,:])
        plot_for_anim(short_data, params, filemarker='sd')
    if args.plot_triangles:
        data, params = get_data(args.input_files[0])
        plot_triangles(data, params)
        plot_simple_conflicts(data, params)
    if args.plot_cpi_conflicts:
        for filename in args.input_files:
            slash = filename.find('/')
            folder = filename[:slash+1]
            if 'Times' in filename:
                time_data, params = get_data(filename)
            elif 'ids' in filename:
                ids_data, cpiparams = read_ids(filename)
            elif 'Conflicts' in filename:
                conflicts_data, cpiparams = read_ids(filename, header_rows=0)
        plot_conflicts(conflicts_data, time_data, ids_data, cpiparams, folder, average=args.average, ax=myax)
    if args.plot_cpi_minority_fraction:
        for filename in args.input_files:
            slash = filename.find('/')
            folder = filename[:slash+1]
            if 'Adj' in filename:
                adj_data, cpiparams = get_data(filename)
            if 'Opns' in filename:
                opns_data, cpiparams = get_data(filename)
            elif 'Times' in filename:
                time_data, params = get_data(filename)
            elif 'ids' in filename:
                ids_data, cpiparams = read_ids(filename)
        plot_timecourse(adj_data, opns_data, time_data, ids_data, cpiparams, folder, cvmp.get_minority_fraction, average=args.average)
    if args.plot_cpi_phase_portrait:
        for filename in args.input_files:
            slash = filename.find('/')
            folder = filename[:slash+1]
            if 'Minorities' in filename:
                minorities_data, cpiparams = read_ids(filename, conversion=float)
            elif 'Conflicts' in filename:
                conflicts_data, cpiparams = read_ids(filename, header_rows=0)
        plot_cpi_phase_portrait(minorities_data, conflicts_data, cpiparams)
    if args.cpi_minorities:
        for filename in args.input_files:
            slash = filename.find('/')
            folder = filename[:slash+1]
            if 'Minorities' in filename:
                minorities_data, cpiparams = get_data(filename)
            elif 'Times' in filename:
                times_data, cpiparams = get_data(filename)
        plot_minorities(minorities_data, times_data, cpiparams, ax=myax)
    if args.plot_single_conflicts:
        for filename in args.input_files:
            if "graphstats" in filename:
                graphdata.append(get_data(filename))
        params = graphdata[0][1]
        plot_single_conflicts(graphdata, params, ax=myax, average=args.average)
    if args.single_minorities:
        for filename in args.input_files:
            if "graphstats" in filename:
                graphdata.append(get_data(filename))
        params = graphdata[0][1]
        plot_single_minorities(graphdata, params, ax=myax, average=args.average)
    if args.dmaps_embeddings:
        for f in args.input_files:
            if 'eigvects' in f:
                eigvects,  params = get_data(f, header_rows=1)
            if 'eigvals' in f:
                eigvals, params = get_data(f, header_rows=1)
            if 'minority' in f:
                label = 'mf'
            else:
                label = 'full'
        vm_dmaps_embedding(eigvals, eigvects, params, label=label)
    if args.cherries:
        ccs = []
        gds = []
        for f in args.input_files:
            if 'cherry' in f:
                cherry_conflicts, g = get_data(f, header_rows=0)
                ccs.append(cherry_conflicts)
            elif 'graphstats' in f:
                graphdata, params = get_data(f, header_rows=1)
                gds.append(graphdata[:,1])
        conflicts_phaseplot(gds, ccs, params)
    if args.squares:
        ccs = []
        gds = []
        for f in args.input_files:
            if 'square' in f:
                cherry_conflicts, g = get_data(f, header_rows=0)
                ccs.append(cherry_conflicts)
            elif 'graphstats' in f:
                graphdata, params = get_data(f, header_rows=1)
                gds.append(graphdata[:,1])
        conflicts_phaseplot(gds, ccs, params)
    if args.timecourses:
        for f in args.input_files:
            if 'cherry' in f:
                cherry_conflicts, g = get_data(f, header_rows=0)
            elif 'graphstats' in f:
                graphdata, params = get_data(f, header_rows=1)
        plot_timecourses(graphdata[:,1], graphdata[:,2], cherry_conflicts, graphdata[:,0], params)

# if args.plot_cpi_phase_portrait:
#         for filename in args.input_files:
#             slash = filename.find('/')
#             folder = filename[:slash+1]
#             if 'Adj' in filename:
#                 adj_data, cpiparams = get_data(filename)
#             if 'Opns' in filename:
#                 opns_data, cpiparams = get_data(filename)
#             elif 'Times' in filename:
#                 time_data, params = get_data(filename)
#             elif 'ids' in filename:
#                 ids_data, cpiparams = read_ids(filename)
#         plot_phase_portrait(adj_data, opns_data, time_data, ids_data, cpiparams, folder, cvmp.get_minority_fraction, cvmp.get_conflicts)
#     myax.legend(fontsize=textsize)
#     plt.show()
