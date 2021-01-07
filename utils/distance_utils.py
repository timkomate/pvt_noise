import math
import numpy as np
import matplotlib.pyplot as plt

def bearing(lat1,lon1,lat2,lon2):
    y = math.sin(lon2-lon1) * math.cos(lat2)
    x = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)* math.cos(lat2) *math.cos(lon2-lon1)
    brng = math.atan2(y,x)
    #math.degrees(brng)
    return brng

def cross_distance(lat1, lon1, lat2, lon2, lat0, lon0, R = 6371.0):
    #d_13 is (angular) distance from start point to third point
    d_13 = angular_distance(
        lat1 = lat1,
        lon1 = lon1, 
        lat2 = lat0,
        lon2 = lon0
    )
    # theta_13 is (initial) bearing from start point to third point
    theta_13 = bearing(
        lat1 = lat1,
        lon1 = lon1,
        lat2 = lat0,
        lon2 = lon0
    )
    # theta_12 is (initial) bearing from start point to end point
    theta_12 = bearing(
        lat1 = lat1,
        lon1 = lon1,
        lat2 = lat2,
        lon2 = lon2
    )
    dXt = math.asin(math.sin(d_13)*math.sin(theta_13-theta_12)) * R
    return dXt

def angular_distance(lat1,lon1,lat2,lon2):
    dlon = lon1 - lon2
    dlat = lat1 - lat2
    a = math.sin(dlat / 2)**2 + math.cos(lat2) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return  c

def spherical_distance(lat1, lon1, lat2, lon2, R=6371.0):
    return angular_distance(lat1,lon1,lat2,lon2) * R

def planar_distance(lat1, lon1, lat2, lon2, lat0, lon0):
    p1 = np.asarray([lat1, lon1])
    p2 = np.asarray([lat2,lon2])
    p3 = np.asarray([lat0, lon0])
    return np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)
