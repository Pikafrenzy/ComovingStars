#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 15:08:13 2025

@author: ambroselo
"""

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
# import os
import pairGraph as pG
from initConditionGenerator import dirCheck
from IDFolder import IDdirCheck
from freeFall import freeFallDifference
from xLimIndex import xLimIndex

# Gala
import gala.dynamics as gd
import gala.potential as gp
# import gala.units as gu

# initial variables
mw = gp.MilkyWayPotential()
dt = 1*u.Myr

def accelEstimate(orbit0,orbit1):
    delta_vel = orbit1.vel-orbit0.vel
    mean_vel = coord.CartesianDifferential(0.5*(orbit1.vel+orbit0.vel).get_d_xyz().to(u.km/u.s))
    delta_pos = orbit1.pos-orbit0.pos
    
    accel = coord.CartesianRepresentation((0.5*delta_vel*mean_vel.norm()/(0.5*delta_pos.norm())).get_d_xyz().to(u.km/u.s**2))
    return accel
    
def trueAccel(orbit):
    accel = mw.acceleration(orbit).to(u.km/u.s**2)
    return accel

def fiducialIntegrate(centre, star0, star1, T):
    wCentre = gd.PhaseSpacePosition(pos = centre.get_Pos(),vel = centre.get_Vel())
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    
    # integrating 
    orbitc = mw.integrate_orbit(wCentre, dt=dt, t1=0, t2 = T)
    orbit0 = mw.integrate_orbit(w0, dt=dt, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt=dt, t1=0, t2 = T)
    
    return orbitc, orbit0, orbit1

def fiducialGraph(ID, T, centre, star0, star1, saveGraphs, dirTime, limitTime):
    freefallLimit = min(limitTime, 100*u.Myr)
    
    orbitc, orbit0, orbit1 = fiducialIntegrate(centre, star0, star1, T)
    star0_freefall_x, star0_freefall_y, star0_freefall_z = freeFallDifference(centre, star0, freefallLimit)
    star1_freefall_x, star1_freefall_y, star1_freefall_z = freeFallDifference(centre, star1, freefallLimit)
    
    aEst = accelEstimate(orbit0, orbit1)
    aTrue = trueAccel(orbitc)
    
    dirCos = pG.getDirCos(orbit0, orbit1)
    
    xMax = int(xLimIndex(orbit0.t[0], limitTime ,dt))

    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'figure.constrained_layout.use':True})
    
    fig = plt.figure(figsize = (12,9))
    
    axDirCos = plt.subplot(4,3,2)
    axDirCos.plot(orbit0.t[0:xMax], (orbit0.t*0)[0:xMax], color = (0.0,0.0,0.0,0.5))
    axDirCos.plot(orbit0.t[0:xMax],dirCos[0:xMax])
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$")
    axDirCos.tick_params('x',labelbottom = False)
    
    axX = plt.subplot(4,3,4)
    axX.plot(orbit0.t[0:xMax],(orbit0.pos.x-orbitc.pos.x)[0:xMax], label = pG.makeLabel(star0.get_Pos(),star0.get_Vel()),color = 'r',lw= 1.0)
    axX.plot(orbit0.t[0:xMax],(orbitc.pos.x-orbitc.pos.x)[0:xMax], label = pG.makeLabel(centre.get_Pos(),centre.get_Vel()),color = 'b',lw= 1.0)
    axX.plot(orbit0.t[0:xMax],(orbit1.pos.x-orbitc.pos.x)[0:xMax], label = pG.makeLabel(star1.get_Pos(),star1.get_Vel()),color = 'g',lw= 1.0)
    axX.plot(orbit0.t[0:len(star0_freefall_x)],star0_freefall_x,label = "Freefall at Red Initial Conditions", color = 'y',linestyle = "--")
    axX.plot(orbit0.t[0:len(star1_freefall_x)],star1_freefall_x,label = "Freefall at Green Initial Conditions",color = 'c',linestyle = "--")
    axX.set_ylabel('Displacement (kpc)')
    axX.set_title("x")
    axX.tick_params('x',labelbottom = False)
    
    fig.legend(loc = 'outside upper left')
    
    axY = plt.subplot(4,3,5, sharey = axX)
    axY.plot(orbit0.t[0:xMax],(orbit0.pos.y-orbitc.pos.y)[0:xMax], color = 'r',lw= 1.0)
    axY.plot(orbit0.t[0:xMax],(orbitc.pos.y-orbitc.pos.y)[0:xMax], color = 'b',lw= 1.0)
    axY.plot(orbit0.t[0:xMax],(orbit1.pos.y-orbitc.pos.y)[0:xMax], color = 'g',lw= 1.0)
    axY.plot(orbit0.t[0:len(star0_freefall_y)],star0_freefall_y,color = 'y',linestyle = "--")
    axY.plot(orbit0.t[0:len(star1_freefall_y)],star1_freefall_y,color = 'c',linestyle = "--")
    axY.set_title("y")
    axY.tick_params('y',labelleft = False)
    axY.tick_params('x',labelbottom = False)
    
    axZ = plt.subplot(4,3,6, sharey = axX)
    axZ.plot(orbit0.t[0:xMax],(orbit0.pos.z-orbitc.pos.z)[0:xMax], color = 'r',lw= 1.0)
    axZ.plot(orbit0.t[0:xMax],(orbitc.pos.z-orbitc.pos.z)[0:xMax], color = 'b',lw= 1.0)
    axZ.plot(orbit0.t[0:xMax],(orbit1.pos.z-orbitc.pos.z)[0:xMax], color = 'g',lw= 1.0)
    axZ.plot(orbit0.t[0:len(star0_freefall_z)],star0_freefall_z,color = 'y',linestyle = "--")
    axZ.plot(orbit0.t[0:len(star1_freefall_z)],star1_freefall_z,color = 'c',linestyle = "--")
    axZ.set_title("z")
    axZ.tick_params('y',labelleft = False)
    axZ.tick_params('x',labelbottom = False)
    
    axVelX = plt.subplot(4,3,7, sharex = axX)
    axVelX.plot(orbit0.t[0:xMax],(orbit0.v_x-orbitc.v_x)[0:xMax], color = 'r',lw= 1.0)
    axVelX.plot(orbit0.t[0:xMax],(orbitc.v_x-orbitc.v_x)[0:xMax], color = 'b',lw= 1.0)
    axVelX.plot(orbit0.t[0:xMax],(orbit1.v_x-orbitc.v_x)[0:xMax], color = 'g',lw= 1.0)
    axVelX.set_title(r"$v_x$")
    axVelX.set_ylabel('Velocity (kpc/Myr)')
    axVelX.tick_params('x',labelbottom = False)

    axVelY = plt.subplot(4,3,8, sharey = axVelX, sharex = axY)
    axVelY.plot(orbit0.t[0:xMax],(orbit0.v_y-orbitc.v_y)[0:xMax], color = 'r',lw= 1.0)
    axVelY.plot(orbit0.t[0:xMax],(orbitc.v_y-orbitc.v_y)[0:xMax], color = 'b',lw= 1.0)
    axVelY.plot(orbit0.t[0:xMax],(orbit1.v_y-orbitc.v_y)[0:xMax], color = 'g',lw= 1.0)
    axVelY.set_title(r"$v_y$")
    axVelY.tick_params('y',labelleft = False)
    axVelY.tick_params('x',labelbottom = False)

    axVelZ = plt.subplot(4,3,9, sharey = axVelX, sharex = axZ)
    axVelZ.plot(orbit0.t[0:xMax],(orbit0.v_z-orbitc.v_z)[0:xMax], color = 'r',lw= 1.0)
    axVelZ.plot(orbit0.t[0:xMax],(orbitc.v_z-orbitc.v_z)[0:xMax], color = 'b',lw= 1.0)
    axVelZ.plot(orbit0.t[0:xMax],(orbit1.v_z-orbitc.v_z)[0:xMax], color = 'g',lw= 1.0)
    axVelZ.set_title(r"$v_z$")
    axVelZ.tick_params('y',labelleft = False)
    axVelZ.tick_params('x',labelbottom = False)
    
    axAccelX = plt.subplot(4,3,10,sharex = axX)
    axAccelX.plot(orbit0.t[0:xMax],aEst[0:xMax].x,color = 'r',label = "Estimated Acceleration")
    axAccelX.plot(orbit0.t[0:xMax],aTrue[0][0:xMax],color = 'g',label = "True Acceleration")
    axAccelX.set_xlabel('Time (Myr)')
    axAccelX.set_title(r"$a_x$")
    axAccelX.set_ylabel('Acceleration($km/s^2$)')
    axAccelX.legend()
    
    axAccelY = plt.subplot(4,3,11,sharex = axY, sharey = axAccelX)
    axAccelY.plot(orbit0.t[0:xMax],aEst[0:xMax].y,color = 'r',label = "Estimated Acceleration")
    axAccelY.plot(orbit0.t[0:xMax],aTrue[1][0:xMax],color = 'g',label = "True Acceleration")
    axAccelY.set_xlabel('Time (Myr)')
    axAccelY.set_title(r"$a_y$")
    axAccelY.tick_params('y',labelleft = False)
    axAccelY.legend()
    
    axAccelZ = plt.subplot(4,3,12,sharex = axZ,sharey = axAccelX)
    axAccelZ.plot(orbit0.t[0:xMax],aEst[0:xMax].z,color = 'r',label = "Estimated Acceleration")
    axAccelZ.plot(orbit0.t[0:xMax],aTrue[2][0:xMax],color = 'g',label = "True Acceleration")
    axAccelZ.set_xlabel('Time (Myr)')
    axAccelZ.set_title(r"$a_z$")
    axAccelZ.tick_params('y',labelleft = False)
    axAccelZ.legend()
    
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