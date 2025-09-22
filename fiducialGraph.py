#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 15:08:13 2025

@author: ambroselo
"""

# import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
# import os
import pairGraph as pG
from initConditionGenerator import dirCheck
from IDFolder import IDdirCheck
from freeFall import freeFallDifference

# Gala
import gala.dynamics as gd
import gala.potential as gp
# import gala.units as gu

# initial variables
mw = gp.MilkyWayPotential()
dt = 1*u.Myr

def fiducialIntegrate(centre, star0, star1, T):
    wCentre = gd.PhaseSpacePosition(pos = centre.get_Pos(),vel = centre.get_Vel())
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    
    # integrating 
    orbitc = mw.integrate_orbit(wCentre, dt=dt, t1=0, t2 = T)
    orbit0 = mw.integrate_orbit(w0, dt=dt, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt=dt, t1=0, t2 = T)
    
    return orbitc, orbit0, orbit1

def fiducialGraph(ID, T, centre, star0, star1, saveGraphs, dirTime):
    orbitc, orbit0, orbit1 = fiducialIntegrate(centre, star0, star1, T)
    star0_freefall_x, star0_freefall_y, star0_freefall_z = freeFallDifference(centre, star0, 100*u.Myr)
    star1_freefall_x, star1_freefall_y, star1_freefall_z = freeFallDifference(centre, star1, 100*u.Myr)
    
    dirCos = pG.getDirCos(orbit0, orbit1)
    
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'figure.constrained_layout.use':True})
    
    fig = plt.figure(figsize = (12,9))
    
    axDirCos = plt.subplot(332)
    axDirCos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDirCos.plot(orbit0.t,dirCos)
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$")
    axDirCos.tick_params('x',labelbottom = False)
    
    axX = plt.subplot(334)
    axX.plot(orbit0.t,orbit0.pos.x-orbitc.pos.x, label = pG.makeLabel(star0.get_Pos(),star0.get_Vel()),color = 'r',lw= 1.0)
    axX.plot(orbit0.t,orbitc.pos.x-orbitc.pos.x, label = pG.makeLabel(centre.get_Pos(),centre.get_Vel()),color = 'b',lw= 1.0)
    axX.plot(orbit0.t,orbit1.pos.x-orbitc.pos.x, label = pG.makeLabel(star1.get_Pos(),star1.get_Vel()),color = 'g',lw= 1.0)
    axX.plot(orbit0.t[0:len(star0_freefall_x)],star0_freefall_x,label = "Freefall at Red Initial Conditions", color = 'y',linestyle = "--")
    axX.plot(orbit0.t[0:len(star1_freefall_x)],star1_freefall_x,label = "Freefall at Green Initial Conditions",color = 'c',linestyle = "--")
    axX.set_ylabel('Displacement (kpc)')
    axX.set_title("x")
    axX.tick_params('x',labelbottom = False)
    
    fig.legend(loc = 'outside upper left')
    
    axY = plt.subplot(335, sharey = axX)
    axY.plot(orbit0.t,orbit0.pos.y-orbitc.pos.y, color = 'r',lw= 1.0)
    axY.plot(orbit0.t,orbitc.pos.y-orbitc.pos.y, color = 'b',lw= 1.0)
    axY.plot(orbit0.t,orbit1.pos.y-orbitc.pos.y, color = 'g',lw= 1.0)
    axY.plot(orbit0.t[0:len(star0_freefall_y)],star0_freefall_y,color = 'y',linestyle = "--")
    axY.plot(orbit0.t[0:len(star1_freefall_y)],star1_freefall_y,color = 'c',linestyle = "--")
    axY.set_title("y")
    axY.tick_params('y',labelleft = False)
    axY.tick_params('x',labelbottom = False)
    
    axZ = plt.subplot(336, sharey = axX)
    axZ.plot(orbit0.t,orbit0.pos.z-orbitc.pos.z, color = 'r',lw= 1.0)
    axZ.plot(orbit0.t,orbitc.pos.z-orbitc.pos.z, color = 'b',lw= 1.0)
    axZ.plot(orbit0.t,orbit1.pos.z-orbitc.pos.z, color = 'g',lw= 1.0)
    axZ.plot(orbit0.t[0:len(star0_freefall_z)],star0_freefall_z,color = 'y',linestyle = "--")
    axZ.plot(orbit0.t[0:len(star1_freefall_z)],star1_freefall_z,color = 'c',linestyle = "--")
    axZ.set_title("z")
    axZ.tick_params('y',labelleft = False)
    axZ.tick_params('x',labelbottom = False)
    
    axVelX = plt.subplot(337)
    axVelX.plot(orbit0.t,orbit0.v_x-orbitc.v_x, color = 'r',lw= 1.0)
    axVelX.plot(orbit0.t,orbitc.v_x-orbitc.v_x, color = 'b',lw= 1.0)
    axVelX.plot(orbit0.t,orbit1.v_x-orbitc.v_x, color = 'g',lw= 1.0)
    axVelX.set_title(r"$v_x$")
    axVelX.set_ylabel('Velocity (kpc/Myr)')
    axVelX.set_xlabel('Time (Myr)')
    
    axVelY = plt.subplot(338, sharey = axVelX)
    axVelY.plot(orbit0.t,orbit0.v_y-orbitc.v_y, color = 'r',lw= 1.0)
    axVelY.plot(orbit0.t,orbitc.v_y-orbitc.v_y, color = 'b',lw= 1.0)
    axVelY.plot(orbit0.t,orbit1.v_y-orbitc.v_y, color = 'g',lw= 1.0)
    axVelY.set_title(r"$v_y$")
    axVelY.tick_params('y',labelleft = False)
    axVelY.set_xlabel('Time (Myr)')

    axVelZ = plt.subplot(339, sharey = axVelX)
    axVelZ.plot(orbit0.t,orbit0.v_z-orbitc.v_z, color = 'r',lw= 1.0)
    axVelZ.plot(orbit0.t,orbitc.v_z-orbitc.v_z, color = 'b',lw= 1.0)
    axVelZ.plot(orbit0.t,orbit1.v_z-orbitc.v_z, color = 'g',lw= 1.0)
    axVelZ.set_title(r"$v_z$")
    axVelZ.tick_params('y',labelleft = False)
    axVelZ.set_xlabel('Time (Myr)')
    
    
    dirPath = dirCheck(saveGraphs, dirTime)
    # saving graphs
    fileName = "Fiducial_and_Star_Pair.png"
    IDdirPath = IDdirCheck(ID,dirPath,saveGraphs)
    pathFiducial = IDdirPath/fileName
    
    if(saveGraphs): 
        plt.savefig(pathFiducial)   
    else:
        plt.show()
    plt.close(fig)