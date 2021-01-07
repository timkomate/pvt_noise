import scipy.io
import numpy as np
import json

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
    
def save_pv(filename,c_branches,distance,model,freqs,gamma,gammaw,
            lat1,lat2,lon1,lon2,nstack,c_zeros = None, f_zeros = None):
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
        "fs": model["freq"].to_numpy(),
        "pvs": model["phase_vel"].to_numpy()
    }
    
    scipy.io.savemat(
        file_name=filename,
        mdict=dd,
        appendmat = True
    )

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