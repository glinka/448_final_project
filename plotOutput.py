import numpy as np
import matplotlib.pyplot as plt
import os

def get_data(fileName, header_rows=1):
    pathToFile = os.path.realpath(fileName)
    f = open(pathToFile, "r")
    params_str = f.readline()
    params = getParams(params_str)
    f.close()
    data = np.genfromtxt(pathToFile, delimiter=",", skip_header=header_rows)
    print params
    return params, data

def get_params(params_str):
    BEGIN = 0
    comma = 1
    #create dict from header, based on key=value format in csv
    params = {}
    while comma > 0:
        comma = params_str.find(",")
        equals = params_str.find("=")
        params[params_str[BEGIN:equals]] = float(params_str[equals+1:comma])
        params_str = params_str[comma+1:]
    params[params_str[BEGIN:equals]] = float(params_str[equals+1:comma])
    #make integer, may not work especially well
    for key in params:
        if(params[key] % 1 == 0):
            params[key] = int(params[key])
    return params

def make_folder(baseName):
    folderCount = 0
    newFolder = baseName + str(folderCount) + '/'
    while os.path.exists(newFolder):
        folderCount = folderCount + 1
        newFolder = baseName + str(folderCount) + '/'
    os.mkdir(newFolder)
    return newFolder

def make_file_name(baseName, params, uniqueID=''):
    fileName = baseName
    for key in params.keys():
        fileName = fileName + '_' + key + '_' + str(params[key])
    if not uniqueID:
        return fileName + '_' + uniqueID
    else
        return fileName

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    parser.add_argument('-cpi', action='store_true', default=False)
    

