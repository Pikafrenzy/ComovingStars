#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 14:56:43 2024

@author: ambroselo
"""

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from Star import Star

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu


rng = np.random.default_rng(137)

stars = []


for i in range(1000):
    
    posR = rng.uniform(0.0,20.0)
    posTheta = rng.uniform(0.0,2*np.pi)
    
    posx = posR*np.cos(posTheta)
    posy = posR*np.sin(posTheta)
    posz = 0.0
    
    vx = 0.0
    vy = rng.normal(220,30)
    vz = 0.0
    
    stars.append(Star(posx,posy,posz,vx,vy,vz))
    