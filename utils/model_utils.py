import pandas as pd
import utils.distance_utils as du
import numpy as np
import warnings

def find_points(lon1,lat1,lon2,lat2):
    lon_mid = np.abs(lon1+lon2)/2
    lat_mid = np.abs(lat1+lat2)/2
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)
    thold2 = du.spherical_distance(lat1,lon1,lat2,lon2) * 0.5
    points = []
    while not points:
        for lon in np.arange(-179.5,179.5+1,1):
            for lat in np.arange(-89.5,89.5+1,1):
                try:
                    dist1 = du.cross_distance(
                        lat1,lon1,lat2,lon2,np.radians(lat),np.radians(lon)
                    )
                    dist2 = du.spherical_distance(
                        np.radians(lat_mid),np.radians(lon_mid),
                        np.radians(lat),np.radians(lon)
                    )
                except:
                    dist1 = np.inf
                    dist2 = np.inf
                
                if(np.abs(dist1) < 60 and dist2 < thold2):
                    #print(lon,lat,dist1, dist2)
                    points.append([lon,lat])
        thold2 *= 1.1
    return np.array(points)

def average_model(points,path,pattern):
    n = points.shape[0]
    pv = []
    for i in np.arange(n):
        lat, lon = points[i,:]
        filename = pattern.format(lon,lat)
        df = pd.read_csv(
            filepath_or_buffer= "{}/{}".format(path,filename),
            delimiter=" ",
            header=None,
            names=["mode","period","phase_vel","group_vel"]
        )
        pv.append(df["phase_vel"].values)

    pv = np.asarray(pv)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        pv = np.nanmean(pv,axis=0)
    data = {
        "period": df["period"].values,
        "freq": 1./df["period"].values,
        "phase_vel": pv
    }
    model = pd.DataFrame(data=data)
    model[model == np.nan] = 2
    model.fillna(2,inplace = True)
    #print(model)
    return model
