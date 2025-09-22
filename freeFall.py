#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 22:32:04 2025

@author: ambroselo
"""

import astropy.units as u
import numpy as np

dt = 1*u.Myr

def freeFallPrediction(star, limit):
    freefall_x = []
    freefall_y = []
    freefall_z = []
    time = 0*u.Myr
    
    while time < limit:
        freefall_x.append((star.get_x()+ time*star.get_Vx()).to_value())
        freefall_y.append((star.get_y()+ time*star.get_Vy()).to_value())
        freefall_z.append((star.get_z()+ time*star.get_Vz()).to_value())
        time = time + dt
    return freefall_x, freefall_y, freefall_z

def freeFallDifference(centre, star, limit):
    centre_freefall_x, centre_freefall_y, centre_freefall_z = freeFallPrediction(centre, limit)
    star_freefall_x, star_freefall_y, star_freefall_z = freeFallPrediction(star, limit)
    difference_freefall_x = np.array(star_freefall_x) - np.array(centre_freefall_x)
    difference_freefall_y = np.array(star_freefall_y) - np.array(centre_freefall_y)
    difference_freefall_z = np.array(star_freefall_z) - np.array(centre_freefall_z)
    return difference_freefall_x, difference_freefall_y, difference_freefall_z