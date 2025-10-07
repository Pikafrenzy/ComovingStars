#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:46:15 2024

@author: ambroselo
"""

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from IDFolder import IDdirCheck
from freeFall import freeFallDifference
from xLimIndex import xLimIndex
# from pathlib import Path
# import os


# Gala
import gala.dynamics as gd
import gala.potential as gp
# import gala.units as gu

# initial variables
mw = gp.MilkyWayPotential()
dt = 1*u.Myr

# truncates a number at 3 decimal places for labels
def trunc3dp(number):
    return np.trunc(number*1e3)/1e3
    
# creates a string from a position and velocity vector of astropy units to print a phase space position as text
# eg. makeLabel([1,2,3]*u.m,[4,5,6]*u.m/u.s) returns
# "Position = [1, 2, 3] m, Velocity = [4, 5, 6] m / s"
def makeLabel(pos, vel):
    label = "Position = "
    label += "["+str(trunc3dp(pos[0].to_value()))+", "+str(trunc3dp(pos[1].to_value()))+", "+str(trunc3dp(pos[2].to_value()))+"] "+pos[0].unit.to_string()
    label += ", Velocity = "
    label += "["+str(trunc3dp(vel[0].to_value()))+", "+str(trunc3dp(vel[1].to_value()))+", "+str(trunc3dp(vel[2].to_value()))+"] "+vel[0].unit.to_string()
    return label

# checks whether the max and min of two arrays are sufficiently near each other for axis limits
def checkSharing(x1, x2, margin):
    diff_min = x2.min()-x1.min()
    diff_max = x2.max()-x1.max()
    share = diff_min.to_value() <= margin and diff_max.to_value() <= margin
    return share

def pairIntegrate(star0, star1, T):
    # creating the PhaseSpacePosition objects from position and velocity vectors
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    
    # integrating orbits
    orbit0 = mw.integrate_orbit(w0, dt = dt, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt = dt, t1=0, t2 = T)
    
    return orbit0, orbit1

def orbitCalcs(orbit0, orbit1):
    diff_pos = orbit1.pos-orbit0.pos
    diff_vel = orbit1.vel-orbit0.vel
    
    # diff_pos.norm()
    
    dirCos = getDirCos(orbit0,orbit1)
    # difference in energy (kinetic, potential, total)
    diffKE = (orbit1.kinetic_energy()-orbit0.kinetic_energy()).to(u.km**2/u.s**2)
    diffPE = (orbit1.potential_energy()-orbit0.potential_energy()).to(u.km**2/u.s**2)
    diffHamiltonian = (orbit1.energy()-orbit0.energy()).to(u.km**2/u.s**2)
    
    return diff_pos, diff_vel, dirCos, diffKE, diffPE, diffHamiltonian

def getDirCos(orbit0, orbit1):
    # difference in each dimension of phase space positions
    diff_pos = orbit1.pos-orbit0.pos

    # calculating the direction cosine
    mean_vel = coord.CartesianDifferential(0.5*(orbit1.vel+orbit0.vel).get_d_xyz().to(u.km/u.s))

    PosDotVel = 0.5*diff_pos.dot(mean_vel)
    dirCos = PosDotVel / (0.5*diff_pos.norm()*mean_vel.norm())
    
    return dirCos

# given two phase space positions
# integrates them using the Milky Way Potential in Gala from Bovy (2015)
# and graphs a series of variables
def pairGraph(ID,star0, star1, T, saveGraphs, dirPath, limitTime):

    orbit0, orbit1 = pairIntegrate(star0, star1, T)
    (diff_pos, diff_vel, dirCos, diffKE,
     diffPE, diffHamiltonian) = orbitCalcs(orbit0, orbit1)
    freefallLimit = min(limitTime, 100*u.Myr)
    
    diff_freefall_x, diff_freefall_y, diff_freefall_z = freeFallDifference(star0, star1, freefallLimit)
    
    xMax = int(xLimIndex(orbit0.t[0], limitTime ,dt))
    
    # creating the initial figure and subfigures
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'figure.constrained_layout.use':True})
    fig = plt.figure(figsize=(16, 12))
    subfigs = fig.subfigures(2,1, height_ratios = [2,4])
    gs0 = gridspec.GridSpec(1,4,figure = subfigs[0])
    gs1 = gridspec.GridSpec(4, 3, figure = subfigs[1])
    
    # filling out the first row of axes
    axY = subfigs[0].add_subplot(gs0[0])
    axZ = subfigs[0].add_subplot(gs0[1],sharey = axY)
    axVelY = subfigs[0].add_subplot(gs0[2])
    axVelZ = subfigs[0].add_subplot(gs0[3],sharey = axVelY)
    
    star0_label = makeLabel(star0.get_Pos(), star0.get_Vel())
    star1_label = makeLabel(star1.get_Pos(), star1.get_Vel())
    
    # plotting the phase space positions of both stars over the integrated period
    axY.plot(orbit0.pos.x[0:xMax],orbit0.pos.y[0:xMax],label=star0_label,linewidth = 0.3, color = 'r')
    axY.plot(orbit1.pos.x[0:xMax],orbit1.pos.y[0:xMax],label=star1_label,linewidth = 0.3, color = 'b')
    axY.set_box_aspect(1)
    axY.axis('equal')
    axY.set_title("Y against X")
    axY.set_ylabel("y (kpc)")
    axY.set_xlabel("x (kpc)")
    
    fig.legend(loc = 'outside upper left')
    
    axZ.plot(orbit0.pos.x[0:xMax],orbit0.pos.z[0:xMax],linewidth = 0.3, color = 'r')
    axZ.plot(orbit1.pos.x[0:xMax],orbit1.pos.z[0:xMax],linewidth = 0.3, color = 'b')
    axZ.set_box_aspect(1)
    axZ.axis('equal')
    axZ.set_ylabel("z (kpc)")
    axZ.set_xlabel("x (kpc)")
    axZ.set_title("Z against X")
    axZ.tick_params('y',labelleft = False)

    axVelY.plot((orbit0.v_x.to(u.km/u.s))[0:xMax],(orbit0.v_y.to(u.km/u.s))[0:xMax],linewidth = 0.3, color = 'r')
    axVelY.plot((orbit1.v_x.to(u.km/u.s))[0:xMax],(orbit1.v_y.to(u.km/u.s))[0:xMax],linewidth = 0.3, color = 'b')
    axVelY.set_box_aspect(1)
    axVelY.axis('equal')
    axVelY.set_title("$v_y$ against $v_x$")
    axVelY.set_ylabel("$v_y$ (km/s)")
    axVelY.set_xlabel("$v_x$ (km/s)")
    
    axVelZ.plot((orbit0.v_x.to(u.km/u.s))[0:xMax],(orbit0.v_z.to(u.km/u.s))[0:xMax],linewidth = 0.3, color = 'r')
    axVelZ.plot((orbit1.v_x.to(u.km/u.s))[0:xMax],(orbit1.v_z.to(u.km/u.s))[0:xMax],linewidth = 0.3, color = 'b')
    axVelZ.set_box_aspect(1)
    axVelZ.axis('equal')
    axVelZ.set_title("$v_z$ against $v_x$")
    axVelZ.set_ylabel("$v_z$ (km/s)")
    axVelZ.set_xlabel("$v_x$ (km/s)")
    axVelZ.tick_params('y',labelleft = False)
    
    # position plots
    axDiffX = subfigs[1].add_subplot(gs1[0,0])
    axDiffX.plot(orbit0.t[0:xMax],diff_pos.x[0:xMax])
    axDiffX.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffX.plot(orbit0.t[0:len(diff_freefall_x)],diff_freefall_x,label = "Freefall", color = 'r')
    axDiffX.set_title("Difference in x")
    axDiffX.set_ylabel("Position (kpc)")
    axDiffX.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_pos.x,diff_pos.y,0.01):
        axDiffY = subfigs[1].add_subplot(gs1[0,1],sharey=axDiffX)
        axDiffY.tick_params('y',labelleft = False)
    else:
        axDiffY = subfigs[1].add_subplot(gs1[0,1])
        axDiffY.set_ylabel("Position (kpc)")
    axDiffY.plot(orbit0.t[0:xMax],diff_pos.y[0:xMax])
    axDiffY.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffY.plot(orbit0.t[0:len(diff_freefall_y)],diff_freefall_y,label = "Freefall", color = 'r')
    axDiffY.set_title("Difference in y")
    axDiffY.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_pos.y,diff_pos.z, 0.01):
        axDiffZ = subfigs[1].add_subplot(gs1[0,2],sharey=axDiffY)
        axDiffZ.tick_params('y',labelleft = False)
    else:
        axDiffZ = fig.add_subplot(gs1[0,2])
        axDiffZ.set_ylabel("Position (kpc)")
    axDiffZ.plot(orbit0.t[0:xMax],diff_pos.z[0:xMax])
    axDiffZ.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffZ.plot(orbit0.t[0:len(diff_freefall_z)],diff_freefall_z,label = "Freefall", color = 'r')
    axDiffZ.set_title("Difference in z")
    axDiffZ.tick_params('x',labelbottom = False)
    
    # velocity plots
    axDiffVx = fig.add_subplot(gs1[1,0],sharex=axDiffX)
    axDiffVx.plot(orbit0.t[0:xMax],diff_vel.d_x[0:xMax])
    axDiffVx.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffVx.set_title("Difference in $v_x$")
    axDiffVx.set_ylabel("Velocity (km/s)")
    axDiffVx.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_vel.d_x,diff_vel.d_y, 0.1):
        axDiffVy = fig.add_subplot(gs1[1,1],sharex = axDiffY, sharey=axDiffVx)
        axDiffVy.tick_params('y',labelleft = False)
    else:
        axDiffVy = fig.add_subplot(gs1[1,1])
        axDiffVy.set_ylabel("Velocity (km/s)")
    axDiffVy.plot(orbit0.t[0:xMax],diff_vel.d_y[0:xMax])
    axDiffVy.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffVy.set_title("Difference in $v_y$")
    axDiffVy.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_vel.d_y,diff_vel.d_z, 0.1):
        axDiffVz = fig.add_subplot(gs1[1,2],sharex = axDiffZ, sharey=axDiffVy)
        axDiffVz.tick_params('y',labelleft = False)
    else:
        axDiffVz = fig.add_subplot(gs1[1,2])
        axDiffVz.set_ylabel("Velocity (km/s)")
    axDiffVz.plot(orbit0.t[0:xMax],diff_vel.d_z[0:xMax])
    axDiffVz.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffVz.set_title("Difference in $v_z$")
    axDiffVz.tick_params('x',labelbottom = False)
    
    # plotting the position displacement magnitude
    axDiffMagPos = fig.add_subplot(gs1[2,0],sharex = axDiffVx)
    axDiffMagPos.plot(orbit0.t[0:xMax],diff_pos.norm()[0:xMax])
    axDiffMagPos.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffMagPos.set_title("Magnitude of difference in position")
    axDiffMagPos.set_ylabel(r"$|\vec x|$ (kpc)")
    axDiffMagPos.tick_params('x',labelbottom = False)
    
    # plotting the velocity displacement magnitude
    axDiffMagVel = fig.add_subplot(gs1[2,2],sharex = axDiffVz)
    axDiffMagVel.plot(orbit0.t[0:xMax], diff_vel.norm()[0:xMax])
    axDiffMagVel.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDiffMagVel.set_title("Magnitude of difference in velocity")
    axDiffMagVel.set_ylabel(r"$|\vec v|$ (km/s)")
    axDiffMagVel.tick_params('x',labelbottom = False)
    
    # plotting the direction cosine
    axDirCos = fig.add_subplot(gs1[2,1],sharex = axDiffVy)
    axDirCos.plot(orbit0.t[0:xMax],dirCos[0:xMax])    
    axDirCos.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$")
    axDirCos.tick_params('x',labelbottom = False)
    
    # energy plots
    axKE = fig.add_subplot(gs1[3,0])
    axKE.plot(orbit0.t[0:xMax],diffKE[0:xMax])
    axKE.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axKE.set_title("Difference in Kinetic Energy")
    axKE.set_ylabel(r"Energy ($\text{km}^2 \text{s}^{-2}$)")
    axKE.set_xlabel("t (Myr)")
    
    axPE = fig.add_subplot(gs1[3,1],sharey = axKE)
    axPE.plot(orbit0.t[0:xMax],diffPE[0:xMax])
    axPE.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axPE.set_title("Difference in Potential Energy")
    axPE.set_xlabel("t (Myr)")
    axPE.tick_params('y',labelleft = False)
    
    axHamiltonian = fig.add_subplot(gs1[3,2])
    axHamiltonian.plot(orbit0.t[0:xMax],diffHamiltonian[0:xMax])
    axHamiltonian.set_title("Difference in Total Energy")
    axHamiltonian.set_xlabel("t (Myr)")

    # saving graphs
    fileName = "Star_Pair.png"
    IDdirPath = IDdirCheck(ID,dirPath,saveGraphs)
    pathPair = IDdirPath/fileName
    
    if(saveGraphs): 
        plt.savefig(pathPair)   
    else:
        plt.show()
    plt.close(fig)
    
    
