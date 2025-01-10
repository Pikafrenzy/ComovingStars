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
import pairGraph

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu


rng = np.random.default_rng(137)

starCentres = []
starPairCount = 12

#toggle for whether graphs are saved automatically to computer
saveGraphs = False

for i in range(starPairCount):
    
    posR = rng.uniform(0.0,20.0)
    posTheta = rng.uniform(0.0,2*np.pi)
    
    posx = posR*np.cos(posTheta)*u.kpc
    posy = posR*np.sin(posTheta)*u.kpc
    posz = 0.0*u.kpc
    
    vx = rng.normal(0.0,30)*u.km/u.s
    vy = rng.normal(220,30)*u.km/u.s
    vz = rng.normal(0.0,30)*u.km/u.s
    
    starCentres.append(Star(posx,posy,posz,vx,vy,vz))

#arrays for the centre of each pair
starsCentreX = []
starsCentreY = []
starsCentreZ = []
starsCentreVx = []
starsCentreVy = []
starsCentreVz = []
for i,star in enumerate(starCentres):
    starsCentreX.append(star.get_x().to_value())
    starsCentreY.append(star.get_y().to_value())
    starsCentreZ.append(star.get_z().to_value())
    starsCentreVx.append(star.get_Vx().to_value())
    starsCentreVy.append(star.get_Vy().to_value())
    starsCentreVz.append(star.get_Vz().to_value())
    
plt.rcParams.update({'font.size': 25})

fig1 = plt.figure(figsize = (20,20),layout = 'constrained')

axCentreXY = plt.subplot(221)
axCentreXY.scatter(starsCentreX,starsCentreY)
axCentreXY.set_title("Y against X of Star Pair Centres")
axCentreXY.set_ylabel("y (kpc)")
axCentreXY.tick_params('x',labelbottom = False)
axCentreXY.tick_params('both',length = 30)

axCentreXZ = plt.subplot(223,sharex=axCentreXY)
axCentreXZ.scatter(starsCentreX,starsCentreZ)
axCentreXZ.set_title("Z against X of Star Pair Centres")
axCentreXZ.set_ylabel("z (kpc)")
axCentreXZ.set_xlabel("x (kpc)")
axCentreXZ.tick_params('both',length = 30)

axCentreVelXY = plt.subplot(222)
axCentreVelXY.scatter(starsCentreVx,starsCentreVy)
axCentreVelXY.set_title("$v_y$ against $v_x$ of Star Pair Centres")
axCentreVelXY.set_ylabel("$v_y$ (km/s)")
axCentreVelXY.tick_params('x',labelbottom = False)
axCentreVelXY.tick_params('both',length = 30)

axCentreVelXZ = plt.subplot(224,sharex=axCentreVelXY)
axCentreVelXZ.scatter(starsCentreVx,starsCentreVz)
axCentreVelXZ.set_title("$v_z$ against $v_x$ of Star Pair Centres")
axCentreVelXZ.set_ylabel("$v_z$ (km/s)")
axCentreVelXZ.set_xlabel("$v_x$ (km/s)")
axCentreVelXZ.tick_params('both',length = 30)

dirPath = "../Plots/"+datetime.now().strftime("%Y%m%d_%H%M%S")
if(saveGraphs): 
   os.mkdir(dirPath)
centreDistributionPath = dirPath+"/StarCentreDistribution_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(centreDistributionPath)
    
    
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
for star in starCentres:
    starPair = pairGenerator(star)
    starPairs.append(starPair)

starX = []
starY = []
starZ = []
starVx = []
starVy = []
starVz = []  

diff_X = []
diff_Y = []
diff_Z = []
diff_Vx = []
diff_Vy = []
diff_Vz = []

for j,starPair in enumerate(starPairs):
    for (i,star) in enumerate(starPair):
        starX.append(starPair[i].get_x().to_value())
        starY.append(starPair[i].get_y().to_value())
        starZ.append(starPair[i].get_z().to_value())
        starVx.append(starPair[i].get_Vx().to_value())
        starVy.append(starPair[i].get_Vy().to_value())
        starVz.append(starPair[i].get_Vz().to_value())
    diff_X.append((starPair[1].get_x()-starPair[0].get_x()).to_value())
    diff_Y.append((starPair[1].get_y()-starPair[0].get_y()).to_value())
    diff_Z.append((starPair[1].get_z()-starPair[0].get_z()).to_value())
    diff_Vx.append((starPair[1].get_Vx()-starPair[0].get_Vx()).to_value())
    diff_Vy.append((starPair[1].get_Vy()-starPair[0].get_Vy()).to_value())
    diff_Vz.append((starPair[1].get_Vz()-starPair[0].get_Vz()).to_value())
    pairGraph.pairGraph(j,starPair[0],starPair[1],saveGraphs,dirPath)
    
plt.rcParams.update({'font.size': 25})

fig2 = plt.figure(figsize = (20,20),layout = 'constrained')

axPairXY = plt.subplot(221)
axPairXY.scatter(diff_X,diff_Y)
axPairXY.set_title("Difference in Y against X")
axPairXY.tick_params('x',labelbottom = False)
axPairXY.set_ylabel("y (kpc)")
axPairXY.tick_params('both',length = 15)

axPairXZ = plt.subplot(223, sharex = axPairXY)
axPairXZ.scatter(diff_X,diff_Z)
axPairXZ.set_title("Difference in Z against X")
axPairXZ.set_ylabel("z (kpc)")
axPairXZ.set_xlabel("x (kpc)")
axPairXZ.tick_params('both',length = 15)

axPairVelXY = plt.subplot(222)
axPairVelXY.scatter(diff_Vx,diff_Vy)
axPairVelXY.set_title("Difference in $v_y$ against $v_x$")
axPairVelXY.tick_params('x',labelbottom = False)
axPairVelXY.set_ylabel("y (km/s)")
axPairVelXY.tick_params('both',length = 15)

axPairVelXZ = plt.subplot(224)
axPairVelXZ.scatter(diff_Vx,diff_Vz)
axPairVelXZ.set_title("Difference in $v_z$ against $v_x$")
axPairVelXZ.set_ylabel("z (km/s)")
axPairVelXZ.set_xlabel("x (km/s)")
axPairVelXZ.tick_params('both',length = 15)

diffPath = dirPath+"/Difference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(diffPath)
  
fig3 = plt.figure(figsize = (20,20),layout = 'constrained')

axPairXY = plt.subplot(221)
axPairXY.scatter(starX,starY)
axPairXY.set_title("Y against X of all Stars")
axPairXY.set_ylabel("y (kpc)")
axPairXY.tick_params('x',labelbottom = False)
axPairXY.tick_params('both',length = 15)

axPairXZ = plt.subplot(223,sharex=axPairXY)
axPairXZ.scatter(starX,starZ)
axPairXZ.set_title("Z against X of all Stars")
axPairXZ.set_ylabel("z (kpc)")
axPairXZ.set_xlabel("x (kpc)")
axPairXZ.tick_params('both',length = 15)

axPairVelXY = plt.subplot(222)
axPairVelXY.scatter(starVx,starVy)
axPairVelXY.set_title("$v_y$ against $v_x$ of all Stars")
axPairVelXY.set_ylabel("$v_y$ (km/s)")
axPairVelXY.tick_params('x',labelbottom = False)
axPairVelXY.tick_params('both',length = 15)

axPairVelXZ = plt.subplot(224,sharex=axPairVelXY)
axPairVelXZ.scatter(starVx,starVz)
axPairVelXZ.set_title("$v_z$ against $v_x$ of all Stars")
axPairVelXZ.set_ylabel("$v_z$ (km/s)")
axPairVelXZ.set_xlabel("$v_x$ (km/s)")
axPairVelXZ.tick_params('both',length = 15)

starDistributionPath = dirPath+"/StarDistribution_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(starDistributionPath)
    