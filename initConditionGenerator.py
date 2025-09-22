#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 14:56:43 2024

@author: ambroselo
"""

# import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
from Star import Star
import pairGraph
from collections import namedtuple
from pathlib import Path

# Gala
import gala.dynamics as gd
import gala.potential as gp
# import gala.units as gu

mw = gp.MilkyWayPotential()

# declaring namedtuple of array of phase-space positions
PosVel = namedtuple('PosVel',['x','y','z','v_x','v_y','v_z'])

# checks if directory already exists; if not, creates it, returns the path
def dirCheck(saveGraphs, dirTime):
    dirPath = Path('Plots/'+dirTime)
    if saveGraphs and not dirPath.is_dir(): 
        os.mkdir(dirPath)
    return dirPath

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
# returns an array of star objects and a namedtuple of type PosVel with 6 arrays inside 
def generateStarCentres(centreCount, rng):
    starCentres = []
    for i in range(centreCount):
        
        posR = rng.uniform(0.0,20.0)
        posTheta = rng.uniform(0.0,2*np.pi)
        
        posx = posR*np.cos(posTheta)
        posy = posR*np.sin(posTheta)
        posz = 0.0
        
        velMean = generateVelocity(posx,posy,posz)
        vx = rng.normal(velMean[0],30)
        vy = rng.normal(velMean[1],30)
        vz = rng.normal(velMean[2],30)
        
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
    
    starCoords = PosVel(starsCentreX, starsCentreY, starsCentreZ, starsCentreVx, starsCentreVy, starsCentreVz)
    return starCentres, starCoords

# graphs the centres generated above
def starCentresGraph(starCoords,saveGraphs,dirTime):
    plt.rcParams.update({'font.size': 25})
    
    # plotting the distribution of the centre of star pairs as four plots in one figure
    fig1 = plt.figure(figsize = (20,20),layout = 'constrained')
    
    axCentreXY = plt.subplot(221)
    axCentreXY.scatter(starCoords.x,starCoords.y)
    axCentreXY.set_title("Y against X of Star Pair Centres")
    axCentreXY.set_ylabel("y (kpc)")
    axCentreXY.tick_params('x',labelbottom = False)
    axCentreXY.tick_params('both',length = 15)
    
    axCentreXZ = plt.subplot(223,sharex=axCentreXY)
    axCentreXZ.scatter(starCoords.x,starCoords.z)
    axCentreXZ.set_title("Z against X of Star Pair Centres")
    axCentreXZ.set_ylabel("z (kpc)")
    axCentreXZ.set_xlabel("x (kpc)")
    axCentreXZ.tick_params('both',length = 15)
    
    axCentreVelXY = plt.subplot(222)
    axCentreVelXY.scatter(starCoords.v_x,starCoords.v_y)
    axCentreVelXY.set_title("$v_y$ against $v_x$ of Star Pair Centres")
    axCentreVelXY.set_ylabel("$v_y$ (km/s)")
    axCentreVelXY.tick_params('x',labelbottom = False)
    axCentreVelXY.tick_params('both',length = 15)
    
    axCentreVelXZ = plt.subplot(224,sharex=axCentreVelXY)
    axCentreVelXZ.scatter(starCoords.v_x,starCoords.v_z)
    axCentreVelXZ.set_title("$v_z$ against $v_x$ of Star Pair Centres")
    axCentreVelXZ.set_ylabel("$v_z$ (km/s)")
    axCentreVelXZ.set_xlabel("$v_x$ (km/s)")
    axCentreVelXZ.tick_params('both',length = 15)
    
    # saves graph
    if(saveGraphs): 
        dirPath = dirCheck(saveGraphs, dirTime)
        fileName = "CentreDistribution_"+dirTime+".png"
        centreDistributionPath = dirPath/fileName
        plt.savefig(centreDistributionPath)

# given a point in phase space
# generates two points equally distanced from the point in all 6 dimensions of phase space
# using a multidimensional Gaussian
def pairGenerator(differenceRNG, starCentre, separationScalar):
    
    xDiff = differenceRNG.normal(0,0.001*separationScalar)
    x1 = starCentre.get_x().to_value() + xDiff
    x2 = starCentre.get_x().to_value() - xDiff
    
    yDiff = differenceRNG.normal(0,0.001*separationScalar)
    y1 = starCentre.get_y().to_value() + yDiff
    y2 = starCentre.get_y().to_value() - yDiff
    
    zDiff = differenceRNG.normal(0,0.001*separationScalar)
    z1 = starCentre.get_z().to_value() + zDiff
    z2 = starCentre.get_z().to_value() - zDiff
    
    vxDiff = differenceRNG.normal(0,0.1*separationScalar)
    vx1 = starCentre.get_Vx().to_value() + vxDiff
    vx2 = starCentre.get_Vx().to_value() - vxDiff
    
    vyDiff = differenceRNG.normal(0,0.1*separationScalar)
    vy1 = starCentre.get_Vy().to_value() + vyDiff
    vy2 = starCentre.get_Vy().to_value() - vyDiff
    
    vzDiff = differenceRNG.normal(0,0.1*separationScalar)
    vz1 = starCentre.get_Vz().to_value() + vzDiff
    vz2 = starCentre.get_Vz().to_value() - vzDiff
    
    star1 = Star(x1,y1,z1,vx1,vy1,vz1)
    star2 = Star(x2,y2,z2,vx2,vy2,vz2)
    return [star1, star2]

# creates all star pairs from the centres given
def starPairsCreate(differenceRNG, starCentres, separationScalar):
    starPairs = []
    for star in starCentres:
        starPair = pairGenerator(differenceRNG, star, separationScalar)
        starPairs.append(starPair)
    return starPairs

def starPairGraphs(starPairs,T, saveGraphs, graphIndividualPairs, graphLimit, dirTime):
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
        if graphIndividualPairs and j < graphLimit:
           dirPath = dirCheck(saveGraphs, dirTime)
           pairGraph.pairGraph(j,starPair[0],starPair[1],T,saveGraphs,dirPath)

    # plotting the displacement between the stars as four plots in one figure 
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

    # saving the plot as a png to the same directory as before
    if(saveGraphs): 
        dirPath = dirCheck(saveGraphs, dirTime)
        fileName = "Difference_"+dirTime+".png"
        diffPath = dirPath/fileName
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
    if(saveGraphs): 
        fileName = "StarDistribution_"+dirTime+".png"
        starDistributionPath = dirPath/fileName
        plt.savefig(starDistributionPath)

# calculating the energy difference between two stars at the starting condition
def energyDifference(star0,star1):
    # converting the Star objects into PhaseSpacePosition objects
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    
    # calculating the energy for each position in the Milky Way Potential
    star0PE = mw.energy(w0)
    star0KE = w0.kinetic_energy()
    star0TotalEnergy = star0PE+star0KE
    
    star1PE = mw.energy(w1)
    star1KE = w1.kinetic_energy()
    star1TotalEnergy = star1PE+star1KE
    
    # I don't know why astropy.Quantity's to_value() returns an array but it does
    return (star1TotalEnergy-star0TotalEnergy).to(u.km**2/u.s**2).to_value()[0]
   
# calculating the difference in Lz between two stars at the starting condition
def LzDifference(star0,star1):
    # converting the Star objects into PhaseSpacePosition objects
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())

    # calculating the angular momentum in z
    Lz0 = w0.angular_momentum()[2]
    Lz1 = w1.angular_momentum()[2]
    return (Lz1-Lz0).to_value(u.kpc**2/u.Myr)
 
def invariantArrays(starPairs):
    energyDifferences = []
    LzDifferences = []
    for starPair in starPairs:
        energyDifferences.append(energyDifference(starPair[0],starPair[1]))
        LzDifferences.append(LzDifference(starPair[0],starPair[1]))
    return energyDifferences, LzDifferences


def invariantGraphs(energyDifferences, LzDifferences, saveGraphs, dirTime):
    plt.rcParams.update({'font.size': 10})
    fig4 = plt.figure(layout = 'constrained')
    
    axDeltaE = plt.subplot(121)
    axDeltaE.hist(energyDifferences,bins = "auto")
    axDeltaE.set_xlim(-500,500)
    axDeltaE.set_title("Initial Energy Difference")
    axDeltaE.set_ylabel("Count")
    axDeltaE.set_xlabel(r"Energy $(\text{km}^2/\text{s}^2)$")
    axDeltaE.plot(energyDifferences[0:12],10*np.arange(1,13),color = 'tab:orange',marker = '.',linestyle = 'None')
    
    axDeltaLz = plt.subplot(122,sharey = axDeltaE)
    axDeltaLz.hist(LzDifferences,bins = "auto")
    axDeltaLz.set_title("Initial Lz Difference")
    axDeltaLz.set_xlabel(r"Angular Momentum $(\text{kpc}^2/\text{Myr})$")
    axDeltaLz.plot(LzDifferences[0:12],10*np.arange(1,13),color = 'tab:orange',marker = '.',linestyle = 'None')
    
    # saving the plot as a png to the same directory as before
    if(saveGraphs): 
        dirPath = dirCheck(saveGraphs, dirTime)
        fileName = "InvariantsDistribution_"+dirTime+".png"
        invariantsDistributionPath = dirPath/fileName
        plt.savefig(invariantsDistributionPath)

