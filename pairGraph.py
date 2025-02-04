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
import os

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu

# initial variables
mw = gp.MilkyWayPotential()
T = 1000*u.Myr

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

# given two phase space positions
# integrates them using the Milky Way Potential in Gala from Bovy (2015)
# and graphs a series of variables
def pairGraph(ID,star0, star1, saveGraphs, dirPath):
    # creating the PhaseSpacePosition objects from position and velocity vectors
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    
    # integrating orbits
    orbit0 = mw.integrate_orbit(w0, dt = 1*u.Myr, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt = 1*u.Myr, t1=0, t2 = T)
    
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
    
    # plotting the phase space positions of both stars over the integrated period
    axY.plot(orbit0.pos.x,orbit0.pos.y,label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axY.plot(orbit1.pos.x,orbit1.pos.y,label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axY.set_box_aspect(1)
    axY.axis('equal')
    axY.set_title("Y against X")
    axY.set_ylabel("y (kpc)")
    axY.set_xlabel("x (kpc)")
    
    fig.legend(loc = 'outside upper left')
    
    axZ.plot(orbit0.pos.x,orbit0.pos.z,label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axZ.plot(orbit1.pos.x,orbit1.pos.z,label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axZ.set_box_aspect(1)
    axZ.axis('equal')
    axZ.set_ylabel("z (kpc)")
    axZ.set_xlabel("x (kpc)")
    axZ.set_title("Z against X")
    axZ.tick_params('y',labelleft = False)

    axVelY.plot(orbit0.v_x.to(u.km/u.s),orbit0.v_y.to(u.km/u.s),label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axVelY.plot(orbit1.v_x.to(u.km/u.s),orbit1.v_y.to(u.km/u.s),label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axVelY.set_box_aspect(1)
    axVelY.axis('equal')
    axVelY.set_title("$v_y$ against $v_x$")
    axVelY.set_ylabel("$v_y$ (km/s)")
    axVelY.set_xlabel("$v_x$ (km/s)")
    
    axVelZ.plot(orbit0.v_x.to(u.km/u.s),orbit0.v_z.to(u.km/u.s),label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axVelZ.plot(orbit1.v_x.to(u.km/u.s),orbit1.v_z.to(u.km/u.s),label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axVelZ.set_box_aspect(1)
    axVelZ.axis('equal')
    axVelZ.set_title("$v_z$ against $v_x$")
    axVelZ.set_ylabel("$v_z$ (km/s)")
    axVelZ.set_xlabel("$v_x$ (km/s)")
    axVelZ.tick_params('y',labelleft = False)
    
    # difference in each dimension of phase space positions
    diff_X = orbit1.pos.x-orbit0.pos.x
    diff_Y = orbit1.pos.y-orbit0.pos.y
    diff_Z = orbit1.pos.z-orbit0.pos.z
    diff_Vx = (orbit1.v_x-orbit0.v_x).to(u.km/u.s)
    diff_Vy = (orbit1.v_y-orbit0.v_y).to(u.km/u.s)
    diff_Vz = (orbit1.v_z-orbit0.v_z).to(u.km/u.s)
    
    # position plots
    axDiffX = subfigs[1].add_subplot(gs1[0,0])
    axDiffX.plot(orbit0.t,diff_X)
    axDiffX.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffX.set_title("Difference in x")
    axDiffX.set_ylabel("Position (kpc)")
    axDiffX.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_X,diff_Y,0.01):
        axDiffY = subfigs[1].add_subplot(gs1[0,1],sharey=axDiffX)
        axDiffY.tick_params('y',labelleft = False)
    else:
        axDiffY = subfigs[1].add_subplot(gs1[0,1])
        axDiffY.set_ylabel("Position (kpc)")
    axDiffY.plot(orbit0.t,diff_Y)
    axDiffY.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffY.set_title("Difference in y")
    axDiffY.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Y,diff_Z, 0.01):
        axDiffZ = subfigs[1].add_subplot(gs1[0,2],sharey=axDiffY)
        axDiffZ.tick_params('y',labelleft = False)
    else:
        axDiffZ = fig.add_subplot(gs1[0,2])
        axDiffZ.set_ylabel("Position (kpc)")
    axDiffZ.plot(orbit0.t,diff_Z)
    axDiffZ.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffZ.set_title("Difference in z")
    axDiffZ.tick_params('x',labelbottom = False)
    
    # velocity plots
    axDiffVx = fig.add_subplot(gs1[1,0],sharex=axDiffX)
    axDiffVx.plot(orbit0.t,diff_Vx)
    axDiffVx.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVx.set_title("Difference in $v_x$")
    axDiffVx.set_ylabel("Velocity (km/s)")
    axDiffVx.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Vx,diff_Vy, 0.1):
        axDiffVy = fig.add_subplot(gs1[1,1],sharex = axDiffY, sharey=axDiffVx)
        axDiffVy.tick_params('y',labelleft = False)
    else:
        axDiffVy = fig.add_subplot(gs1[1,1])
        axDiffVy.set_ylabel("Velocity (km/s)")
    axDiffVy.plot(orbit0.t,diff_Vy)
    axDiffVy.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVy.set_title("Difference in $v_y$")
    axDiffVy.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Vy,diff_Vz, 0.1):
        axDiffVz = fig.add_subplot(gs1[1,2],sharex = axDiffZ, sharey=axDiffVy)
        axDiffVz.tick_params('y',labelleft = False)
    else:
        axDiffVz = fig.add_subplot(gs1[1,2])
        axDiffVz.set_ylabel("Velocity (km/s)")
    axDiffVz.plot(orbit0.t,diff_Vz)
    axDiffVz.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVz.set_title("Difference in $v_z$")
    axDiffVz.tick_params('x',labelbottom = False)
    
    # magnitude of the position displacement between the stars
    magPosDifference = np.sqrt(diff_X**2+diff_Y**2+diff_Z**2)
    
    # plotting the position displacement magnitude
    axDiffMagPos = fig.add_subplot(gs1[2,0],sharex = axDiffVx)
    axDiffMagPos.plot(orbit0.t,magPosDifference)
    axDiffMagPos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagPos.set_title("Magnitude of difference in position")
    axDiffMagPos.set_ylabel(r"$|\vec x|$ (kpc)")
    axDiffMagPos.tick_params('x',labelbottom = False)
    
    # magnitude of the velocity displacement between the stars
    magVelDifference = np.sqrt(diff_Vx**2+diff_Vy**2+diff_Vz**2)
    
    # plotting the velocity displacement magnitude
    axDiffMagVel = fig.add_subplot(gs1[2,2],sharex = axDiffVz)
    axDiffMagVel.plot(orbit0.t,magVelDifference)
    axDiffMagVel.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagVel.set_title("Magnitude of difference in velocity")
    axDiffMagVel.set_ylabel(r"$|\vec v|$ (km/s)")
    axDiffMagVel.tick_params('x',labelbottom = False)
    
    # calculating the direction cosine
    mean_Vx = 0.5*(orbit1.v_x+orbit0.v_x).to(u.km/u.s)
    mean_Vy = 0.5*(orbit1.v_y+orbit0.v_y).to(u.km/u.s)
    mean_Vz = 0.5*(orbit1.v_z+orbit0.v_z).to(u.km/u.s)
    
    magVelMean = np.sqrt(mean_Vx**2+mean_Vy**2+mean_Vz**2)
    
    dirCos = []
    for (i, pos) in enumerate(diff_X):
        PosDotVel = 0.5*diff_X[i]*mean_Vx[i]+0.5*diff_Y[i]*mean_Vy[i]+0.5*diff_Z[i]*mean_Vz[i]
        dirCos.append(PosDotVel / (0.5*magPosDifference[i] * magVelMean[i]))
    
    # plotting the direction cosine
    axDirCos = fig.add_subplot(gs1[2,1],sharex = axDiffVy)
    axDirCos.plot(orbit0.t,dirCos)    
    axDirCos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$")
    axDirCos.tick_params('x',labelbottom = False)
    
    # difference in energy (kinetic, potential, total)
    diffKE = (orbit1.kinetic_energy()-orbit0.kinetic_energy()).to(u.km**2/u.s**2)
    diffPE = (orbit1.potential_energy()-orbit0.potential_energy()).to(u.km**2/u.s**2)
    diffHamiltonian = (orbit1.energy()-orbit0.energy()).to(u.km**2/u.s**2)
    
    # energy plots
    axKE = fig.add_subplot(gs1[3,0])
    axKE.plot(orbit0.t,diffKE)
    axKE.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axKE.set_title("Difference in Kinetic Energy")
    axKE.set_ylabel(r"Energy ($\text{km}^2 \text{s}^{-2}$)")
    axKE.set_xlabel("t (Myr)")
    
    axPE = fig.add_subplot(gs1[3,1],sharey = axKE)
    axPE.plot(orbit0.t,diffPE)
    axPE.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axPE.set_title("Difference in Potential Energy")
    axPE.set_xlabel("t (Myr)")
    axPE.tick_params('y',labelleft = False)
    
    axHamiltonian = fig.add_subplot(gs1[3,2])
    axHamiltonian.plot(orbit0.t,diffHamiltonian)
    axHamiltonian.set_title("Difference in Total Energy")
    axHamiltonian.set_xlabel("t (Myr)")

    # saving graphs
    pathPair = dirPath+"/Star_Pair_"+str(ID)+".png"
    
    if(saveGraphs): 
        plt.savefig(pathPair)   
    else:
        plt.show()
    plt.close(fig)
    
    
    
