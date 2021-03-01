import pandas as pd
import utils.distance_utils as du
import numpy as np
import warnings
import matplotlib.pyplot as plt

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
        pv = np.nanmedian(pv,axis=0)

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

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]