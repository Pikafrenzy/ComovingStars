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
saveGraphs = False

for i in range(1000):
    
    posR = rng.uniform(0.0,20.0)
    posTheta = rng.uniform(0.0,2*np.pi)
    
    posx = posR*np.cos(posTheta)*u.kpc
    posy = posR*np.sin(posTheta)*u.kpc
    posz = 0.0*u.kpc
    
    vx = rng.normal(0.0,30)*u.km/u.s
    vy = rng.normal(220,30)*u.km/u.s
    vz = rng.normal(0.0,30)*u.km/u.s
    
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
axXY.tick_params('x',labelbottom = False)
axXY.tick_params('both',length = 30)

axXZ = plt.subplot(223,sharex=axXY)
axXZ.scatter(starsX,starsZ)
axXZ.set_title("Z against X")
axXZ.set_ylabel("z (kpc)")
axXZ.set_xlabel("x (kpc)")
axXZ.tick_params('both',length = 30)

axVelXY = plt.subplot(222)
axVelXY.scatter(starsVx,starsVy)
axVelXY.set_title("$v_y$ against $v_x$")
axVelXY.set_ylabel("$v_y$ (km/s)")
axVelXY.tick_params('x',labelbottom = False)
axVelXY.tick_params('both',length = 30)

axVelXZ = plt.subplot(224,sharex=axVelXY)
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
    
    
def pairGenerator(starCentre):
    
    xDiff = rng.normal(0,0.001)
    x1 = starCentre.get_x().to_value() + xDiff
    x2 = starCentre.get_x().to_value() - xDiff
    
    yDiff = rng.normal(0,0.001)
    y1 = starCentre.get_y().to_value() + yDiff
    y2 = starCentre.get_y().to_value() - yDiff
    
    zDiff = rng.normal(0,0.001)
    z1 = starCentre.get_z().to_value() + zDiff
    z2 = starCentre.get_z().to_value() - zDiff
    
    vxDiff = rng.normal(0,0.1)
    vx1 = starCentre.get_Vx().to_value() + vxDiff
    vx2 = starCentre.get_Vx().to_value() - vxDiff
    
    vyDiff = rng.normal(0,0.1)
    vy1 = starCentre.get_Vy().to_value() + vyDiff
    vy2 = starCentre.get_Vy().to_value() - vyDiff
    
    vzDiff = rng.normal(0,0.1)
    vz1 = starCentre.get_Vz().to_value() + vzDiff
    vz2 = starCentre.get_Vz().to_value() - vzDiff
    
    star1 = Star(x1,y1,z1,vx1,vy1,vz1)
    star2 = Star(x2,y2,z2,vx2,vy2,vz2)
    return [star1, star2]

starPairs = []
for star in stars:
    starPair = pairGenerator(star)
    starPairs.append(starPair)

diff_X = []
diff_Y = []
diff_Z = []
diff_Vx = []
diff_Vy = []
diff_Vz = []

for starPair in starPairs:
    diff_X.append((starPair[1].get_x()-starPair[0].get_x()).to_value())
    diff_Y.append((starPair[1].get_y()-starPair[0].get_y()).to_value())
    diff_Z.append((starPair[1].get_z()-starPair[0].get_z()).to_value())
    diff_Vx.append((starPair[1].get_Vx()-starPair[0].get_Vx()).to_value())
    diff_Vy.append((starPair[1].get_Vy()-starPair[0].get_Vy()).to_value())
    diff_Vz.append((starPair[1].get_Vz()-starPair[0].get_Vz()).to_value())
    
fig2 = plt.figure(figsize = (40,40),layout = 'constrained')

axPairXY = plt.subplot(221)
axPairXY.scatter(diff_X,diff_Y)
axPairXY.set_title("Difference in Y against X")
axPairXY.tick_params('x',labelbottom = False)
axPairXY.set_ylabel("y (kpc)")
axPairXY.tick_params('both',length = 30)

axPairXZ = plt.subplot(223, sharex = axPairXY)
axPairXZ.scatter(diff_X,diff_Z)
axPairXZ.set_title("Difference in Z against X")
axPairXZ.set_ylabel("z (kpc)")
axPairXZ.set_xlabel("x (kpc)")
axPairXZ.tick_params('both',length = 30)

axPairVelXY = plt.subplot(222)
axPairVelXY.scatter(diff_Vx,diff_Vy)
axPairVelXY.set_title("Difference in $v_y$ against $v_x$")
axPairVelXY.tick_params('x',labelbottom = False)
axPairVelXY.set_ylabel("y (km/s)")
axPairVelXY.tick_params('both',length = 30)

axPairVelXZ = plt.subplot(224)
axPairVelXZ.scatter(diff_Vx,diff_Vz)
axPairVelXZ.set_title("Difference in $v_z$ against $v_x$")
axPairVelXZ.set_ylabel("z (km/s)")
axPairVelXZ.set_xlabel("x (km/s)")
axPairVelXZ.tick_params('both',length = 30)

diffPath = dirPath+"/Difference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(diffPath)
