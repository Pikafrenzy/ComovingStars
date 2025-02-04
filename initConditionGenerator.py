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
import time
from Star import Star
import pairGraph

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu

# notes the start time of running the program to determine how long it ran for
startTime = time.time()

# a seeded RNG to ensure consistent data when running multiple times
# TODO: separate out rng to several different rng generators
# such that increasing star pair count doesn't give different star pairs
rng = np.random.default_rng(137) 


starCentres = []
starPairCount = 12

# toggle for whether graphs are saved automatically to computer
saveGraphs = False

# taking a position,
# calculates the cross product of the unit vector in the z axis the normalised position vector
# and returns it with a magnitude of 220

def generateVelocity(posx,posy,posz):
    magPos = np.sqrt(posx**2+posy**2+posz**2)
    normPosX = posx/magPos
    normPosY = posy/magPos
    normPosZ = posz/magPos
    normPosVec = [normPosX,normPosY,normPosZ]
    zHat = [0,0,1]
    vel = np.cross(zHat,normPosVec)*220
    return vel

# generates an array of points in phase space which represent the centres of a pair of stars
for i in range(starPairCount):
    
    posR = rng.uniform(0.0,20.0)
    posTheta = rng.uniform(0.0,2*np.pi)
    
    posx = posR*np.cos(posTheta)*u.kpc
    posy = posR*np.sin(posTheta)*u.kpc
    posz = 0.0*u.kpc
    
    velMean = generateVelocity(posx.to_value(),posy.to_value(),posz.to_value())
    vx = rng.normal(velMean[0],30)*u.km/u.s
    vy = rng.normal(velMean[1],30)*u.km/u.s
    vz = rng.normal(velMean[2],30)*u.km/u.s
    
    starCentres.append(Star(posx,posy,posz,vx,vy,vz))

# creating and populating arrays containing the centre of each pair
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

# plotting the distribution of the centre of star pairs as four plots in one figure
fig1 = plt.figure(figsize = (20,20),layout = 'constrained')

axCentreXY = plt.subplot(221)
axCentreXY.scatter(starsCentreX,starsCentreY)
axCentreXY.set_title("Y against X of Star Pair Centres")
axCentreXY.set_ylabel("y (kpc)")
axCentreXY.tick_params('x',labelbottom = False)
axCentreXY.tick_params('both',length = 15)

axCentreXZ = plt.subplot(223,sharex=axCentreXY)
axCentreXZ.scatter(starsCentreX,starsCentreZ)
axCentreXZ.set_title("Z against X of Star Pair Centres")
axCentreXZ.set_ylabel("z (kpc)")
axCentreXZ.set_xlabel("x (kpc)")
axCentreXZ.tick_params('both',length = 15)

axCentreVelXY = plt.subplot(222)
axCentreVelXY.scatter(starsCentreVx,starsCentreVy)
axCentreVelXY.set_title("$v_y$ against $v_x$ of Star Pair Centres")
axCentreVelXY.set_ylabel("$v_y$ (km/s)")
axCentreVelXY.tick_params('x',labelbottom = False)
axCentreVelXY.tick_params('both',length = 15)

axCentreVelXZ = plt.subplot(224,sharex=axCentreVelXY)
axCentreVelXZ.scatter(starsCentreVx,starsCentreVz)
axCentreVelXZ.set_title("$v_z$ against $v_x$ of Star Pair Centres")
axCentreVelXZ.set_ylabel("$v_z$ (km/s)")
axCentreVelXZ.set_xlabel("$v_x$ (km/s)")
axCentreVelXZ.tick_params('both',length = 15)

# creates the directory for all the plots, and saves the above plot to it as a png file
dirPath = "../Plots/"+datetime.now().strftime("%Y%m%d_%H%M%S")
if(saveGraphs): 
   os.mkdir(dirPath)
centreDistributionPath = dirPath+"/CentreDistribution_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(centreDistributionPath)
    
    
# given a point in phase space
# generates two points equally distanced from the point in all 6 dimensions of phase space
# using a multidimensional Gaussian
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

# creating and populating arrays containing the locations of each pair and their differences in phase space
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

for (j,starPair) in enumerate(starPairs):
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
    # graphing each pair of stars
    pairGraph.pairGraph(j,starPair[0],starPair[1],saveGraphs,dirPath)
    
    
# plotting the displacement between the stars as four plots in one figure 

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

# saving the plot as a png to the same directory as before
diffPath = dirPath+"/Difference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(diffPath)
  
# plotting the distribution of the stars as four plots in one figure     

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

# saving the plot as a png to the same directory as before
starDistributionPath = dirPath+"/StarDistribution_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(starDistributionPath)
    
# outputs the time it took the program to run
endTime = time.time()
print("This took " + str(endTime-startTime) + " seconds")