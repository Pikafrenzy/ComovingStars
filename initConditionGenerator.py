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
import os
from datetime import datetime
from Star import Star

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu


rng = np.random.default_rng(137)

stars = []

#toggle for whether graphs are saved automatically to computer
saveGraphs = True

for i in range(1000):
    
    posR = rng.uniform(0.0,20.0)
    posTheta = rng.uniform(0.0,2*np.pi)
    
    posx = posR*np.cos(posTheta)*u.kpc
    posy = posR*np.sin(posTheta)*u.kpc
    posz = 0.0*u.kpc
    
    vx = 0.0*u.km/u.s
    vy = rng.normal(220,30)*u.km/u.s
    vz = 0.0*u.km/u.s
    
    stars.append(Star(posx,posy,posz,vx,vy,vz))

starsX = []
starsY = []
starsZ = []
starsVx = []
starsVy = []
starsVz = []
for i,star in enumerate(stars):
    starsX.append(star.get_x().to_value())
    starsY.append(star.get_y().to_value())
    starsZ.append(star.get_z().to_value())
    starsVx.append(star.get_Vx().to_value())
    starsVy.append(star.get_Vy().to_value())
    starsVz.append(star.get_Vz().to_value())
    
plt.rcParams.update({'font.size': 50})

fig1 = plt.figure(figsize = (40,40),layout = 'constrained')

axXY = plt.subplot(221)
axXY.scatter(starsX,starsY)
axXY.set_title("Y against X")
axXY.set_ylabel("y (kpc)")
axXY.set_xlabel("x (kpc)")
axXY.tick_params('both',length = 30)

axXZ = plt.subplot(222)
axXZ.scatter(starsX,starsZ)
axXZ.set_title("Z against X")
axXZ.set_ylabel("z (kpc)")
axXZ.set_xlabel("x (kpc)")
axXZ.tick_params('both',length = 30)

axVelXY = plt.subplot(223)
axVelXY.scatter(starsVx,starsVy)
axVelXY.set_title("$v_y$ against $v_x$")
axVelXY.set_ylabel("$v_y$ (km/s)")
axVelXY.set_xlabel("$v_x$ (km/s)")
axVelXY.tick_params('both',length = 30)

axVelXZ = plt.subplot(224)
axVelXZ.scatter(starsVx,starsVz)
axVelXZ.set_title("$v_z$ against $v_x$")
axVelXZ.set_ylabel("$v_z$ (km/s)")
axVelXZ.set_xlabel("$v_x$ (km/s)")
axVelXZ.tick_params('both',length = 30)

dirPath = "../Plots/"+datetime.now().strftime("%Y%m%d_%H%M%S")
if(saveGraphs): 
   os.mkdir(dirPath)
distPath = dirPath+"/Distribution_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(distPath)