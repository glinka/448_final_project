import calcVotingModelProps as cvmp
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

def read_ids(filename, header_rows=1, **kwargs):
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
            ids[-1].append(int(myline[0:comma_loc]))
            myline = myline[comma_loc+1:]
    return np.array(ids), params


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
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if average:
        ax.plot(time_data, final_ydata)
    else:
        ax.hold(True)
        for i in range(starting_nvms):
            npoints = len(final_ydata[i])
            ax.plot(time_data[:npoints], final_ydata[i])
    plt.show()

def plot_conflicts(conflicts_data, time_data, ids_data, params, folder, average=False, ax=''):
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
            final_ydata.append(np.average([ydata[str(j)][i] for j in ids_data[i]]))
    else:
        for key in ydata.keys():
            final_ydata.append(ydata[key])
    final_ydata = np.array(final_ydata)
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if average:
        ax.plot(time_data, final_ydata)
    else:
        ax.hold(True)
        for i in range(starting_nvms):
            npoints = len(final_ydata[i])
            ax.plot(time_data[:npoints], final_ydata[i])
    plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('-cpic', '--plot-cpi-conflicts', action='store_true', default=False)
    parser.add_argument('-cpimf', '--plot-cpi-minority-fraction', action='store_true', default=False)
    parser.add_argument('-cpipp', '--plot-cpi-phase-portrait', action='store_true', default=False)
    parser.add_argument('-avg', '--average', action='store_true', default=False)
    args = parser.parse_args()
    #should have one of each file: CPIAdj_... and CPIOpns_...
    #which will have the same header (and thus same params)
    for filename in args.input_files:
        slash = filename.find('/')
        folder = filename[:slash+1]
        # if 'Adj' in filename:
        #     adj_data, params = get_data(filename)
        if 'Opns' in filename:
            opns_data, params = get_data(filename)
        elif 'Times' in filename:
            time_data, params = get_data(filename)
        elif 'ids' in filename:
            ids_data, params = read_ids(filename)
        elif 'Conflicts' in filename:
            conflicts_data, params = read_ids(filename, header_rows=0)
    # if args.plot_cpi_conflicts:
    #     plot_timecourse(adj_data, opns_data, time_data, ids_data, params, folder, cvmp.get_conflicts, average=args.average)
    if args.plot_cpi_conflicts:
        plot_conflicts(conflicts_data, time_data, ids_data, params, folder, average=args.average)
    if args.plot_cpi_minority_fraction:
        plot_timecourse(adj_data, opns_data, time_data, ids_data, params, folder, cvmp.get_minority_fraction, average=args.averager)
    if args.plot_cpi_phase_portrait:
        plot_phase_portrait(adj_data, opns_data, time_data, ids_data, params, folder, cvmp.get_minority_fraction, cvmp.get_conflicts)
