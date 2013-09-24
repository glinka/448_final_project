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
    starting_nvms = params['nVms']
    n = params['n']
    n_data = time_data.shape[0]
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ydata = []
    for i in range(n_data):
        n_vms = int(nvms_data[i] + 0.5)
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
            ax.plot(time_data, ydata[:,k])
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
        if 'Adj' in filename:
            adj_data, params = get_data(filename)
        elif 'Opns' in filename:
            opns_data, params = get_data(filename)
        elif 'Times' in filename:
            time_data, params = get_data(filename)
        elif 'Nvms' in filename:
            nvms_data, params = get_data(filename)
    if args.plot_cpi_conflicts:
        plot_timecourse(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_conflicts, average=args.average)
    if args.plot_cpi_minority_fraction:
        plot_timecourse(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_minority_fraction, average=args.average)
    if args.plot_cpi_phase_portrait:
        plot_phase_portrait(adj_data, opns_data, time_data, nvms_data, params, cvmp.get_minority_fraction, cvmp.get_conflicts)
