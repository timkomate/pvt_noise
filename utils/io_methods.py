import scipy.io
import numpy as np
import pandas as pd
import json
import os

def read_json(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return data

def read_measured_data(path):
    matfile = scipy.io.loadmat(
        file_name = path,
        squeeze_me = True
    )
    distance = matfile["Dist"] / 1000 #m to km
    dt = matfile["dtnew"]
    t = matfile["lagsx1x2"]
    ccf = matfile["cross12"]
    nstack = matfile["nstack"]
    s1 = matfile["Station1"]
    s2 = matfile["Station2"]
    try:
        cutvec = matfile["cutvec"]
        return [t,ccf,distance,dt,nstack,s1,s2,cutvec]
    except:
        n1 = matfile["Network1"]
        n2 = matfile["Network2"]
        return [t,ccf,distance,dt,nstack,n1,s1,n2,s2]

def save_bg_model(fname,model):
    model.to_csv(
        path_or_buf = fname,
        sep = " ",
        columns = ["period", "phase_vel", "phase_vel"],
        header = False,
        index = False,
        float_format="%.2f"
    )
    
def save_pv(filename,c_branches,distance,model,freqs,gamma,gammaw,
            lat1,lat2,lon1,lon2,nstack, n1, s1, n2, s2,
            c_zeros = None, f_zeros = None):
    dd = {
        "c_zeros": c_zeros,
        "f_zeros": f_zeros,
        "crayan": c_branches,
        "Dist": distance,
        "frayan": freqs,
        "gamma1": gamma,
        "gammaw": gammaw,
        "lat1": lat1,
        "lat2": lat2,
        "lon1": lon1,
        "lon2": lon2,
        "nstack": nstack,
        "Network1": n1,
        "Station1": s1,
        "Network2": n2,
        "Station2": s2,
        "fs": model["freq"].to_numpy(),
        "pvs": model["phase_vel"].to_numpy()
    }
    
    scipy.io.savemat(
        file_name=filename,
        mdict=dd,
        appendmat = True
    )

def save_pv_format(save_path,filename,c_branches,distance,model,freqs,
            gamma,gammaw,lat1,lat2,lon1,lon2,nstack, n1, s1, n2, s2,
            c_zeros = None, f_zeros = None):
    dd = {
        "c_zeros": c_zeros,
        "f_zeros": f_zeros,
        "crayan": c_branches,
        "Dist": distance,
        "frayan": freqs,
        "gamma1": gamma,
        "gammaw": gammaw,
        "lat1": lat1,
        "lat2": lat2,
        "lon1": lon1,
        "lon2": lon2,
        "nstack": nstack,
        "Network1": n1,
        "Station1": s1,
        "Network2": n2,
        "Station2": s2,
        "fs": model["freq"].to_numpy(),
        "pvs": model["phase_vel"].to_numpy()
    }

    folder = "{}/{}-{}/".format(save_path,s1,s2)
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    scipy.io.savemat(
        file_name="{}/{}".format(folder,filename),
        mdict=dd,
        appendmat = True
    )
    print("{}/{} saved".format(folder,filename))

def save_results_ascii(freqs,pv,filename):
    shape = pv.shape
    output = open(filename,"w")
    for i in np.arange(shape[0]):
        for j in np.arange(shape[1]):
            output.write("{} {}\n".format(freqs[j], pv[i,j]))
        output.write(">\n")
    output.close()

def save_synthetic_ccf(t,ccf,dt,distance):
    dd = {
        "Station1": "S1",
        "Station2": "S2",
        "Dist": distance,
        "dtnew": dt,
        "nstack": 1,
        "lagsx1x2": t,
        "cross12": ccf
    }
    
    scipy.io.savemat(
        "./synthetic.mat",
        dd
    )