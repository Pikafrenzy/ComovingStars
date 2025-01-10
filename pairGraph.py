#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:46:15 2024

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

#initial variables
mw = gp.MilkyWayPotential()
T = 500*u.Myr
def trunc3dp(number):
    return np.trunc(number*1e3)/1e3
    
def makeLabel(vec3):
    label = "["+str(trunc3dp(vec3[0].to_value()))+", "+str(trunc3dp(vec3[1].to_value()))+", "+str(trunc3dp(vec3[2].to_value()))+"] "+vec3[0].unit.to_string()
    return label

def pairGraph(ID,star0, star1, saveGraphs, dirPath):
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    orbit0 = mw.integrate_orbit(w0, dt = 1*u.Myr, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt = 1*u.Myr, t1=0, t2 = T)
    
    plt.rcParams.update({'font.size': 10})
    
    fig1 = plt.figure(layout = "constrained")
    ax1Y = plt.subplot(121)
    ax1Z = plt.subplot(122,sharey = ax1Y)
    
    # fig1, ax1 = plt.subplots(1,2, layout="constrained",sharey = True)
    ax1Y.plot(orbit0.pos.x,orbit0.pos.y,label="Position = "+makeLabel(star0.get_Pos()),linewidth = 0.3, color = 'r')
    ax1Y.plot(orbit1.pos.x,orbit1.pos.y,label="Position = "+makeLabel(star1.get_Pos()),linewidth = 0.3, color = 'b')
    ax1Y.legend(fontsize=10, bbox_to_anchor=(0.9,-0.2))
    ax1Y.set_box_aspect(1)
    ax1Y.set_title("Y against X")
    # ax1Y.set_box_aspect(1)
    
    ax1Z.plot(orbit0.pos.x,orbit0.pos.z,label="Position = "+makeLabel(star0.get_Pos()),linewidth = 0.3, color = 'r')
    ax1Z.plot(orbit1.pos.x,orbit1.pos.z,label="Position = "+makeLabel(star1.get_Pos()),linewidth = 0.3, color = 'b')
    ax1Z.set_box_aspect(1)
    ax1Z.set_title("Z against X")
    ax1Z.tick_params('y',labelleft = False)
    # ax1Z.set_box_aspect(1)
    
    dirPathPair = dirPath+"/"+str(ID)
    posPath = dirPathPair+"/Position_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
    if(saveGraphs): 
        os.mkdir(dirPathPair)
        plt.savefig(posPath)
    
    plt.rcParams.update({'font.size': 28})
    
    diff_X = orbit1.pos.x-orbit0.pos.x
    diff_Y = orbit1.pos.y-orbit0.pos.y
    diff_Z = orbit1.pos.z-orbit0.pos.z
    diff_Vx = (orbit1.v_x-orbit0.v_x).to(u.km/u.s)
    diff_Vy = (orbit1.v_y-orbit0.v_y).to(u.km/u.s)
    diff_Vz = (orbit1.v_z-orbit0.v_z).to(u.km/u.s)
    
    fig2 = plt.figure(figsize = (30,10), layout = "constrained")
    
    #position plots
    axDiffX = plt.subplot(331)
    axDiffX.plot(orbit0.t,diff_X)
    axDiffX.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffX.set_title("Difference in x")
    axDiffX.set_ylabel("Position (kpc)")
    axDiffX.tick_params('x',labelbottom = False)
    axDiffX.tick_params('both',length = 12)
    
    def checkSharing(x1, x2, margin):
        diff_min = x2.min()-x1.min()
        diff_max = x2.max()-x1.max()
        share = diff_min.to_value() <= margin and diff_max.to_value() <= margin
        return share
    
    if checkSharing(diff_X,diff_Y,0.01):
        axDiffY = plt.subplot(332,sharey=axDiffX)
        axDiffY.tick_params('y',labelleft = False)
    else:
        axDiffY = plt.subplot(332)
        axDiffY.set_ylabel("Position (kpc)")
    axDiffY.plot(orbit0.t,diff_Y)
    axDiffY.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffY.set_title("Difference in y")
    axDiffY.tick_params('x',labelbottom = False)
    axDiffY.tick_params('both',length = 12)
    
    if checkSharing(diff_Y,diff_Z, 0.01):
        axDiffZ = plt.subplot(333,sharey=axDiffY)
        axDiffZ.tick_params('y',labelleft = False)
    else:
        axDiffZ = plt.subplot(333)
        axDiffZ.set_ylabel("Position (kpc)")
    axDiffZ.plot(orbit0.t,diff_Z)
    axDiffZ.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffZ.set_title("Difference in z")
    axDiffZ.tick_params('x',labelbottom = False)
    axDiffZ.tick_params('both',length = 12)
    
    #velocity plots
    axDiffVx = plt.subplot(334,sharex=axDiffX)
    axDiffVx.plot(orbit0.t,diff_Vx)
    axDiffVx.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVx.set_title("Difference in $v_x$")
    axDiffVx.set_ylabel("Velocity (km/s)")
    axDiffVx.tick_params('x',labelbottom = False)
    axDiffVx.tick_params('both',length = 12)
    
    if checkSharing(diff_Vx,diff_Vy, 0.1):
        axDiffVy = plt.subplot(335,sharex = axDiffY, sharey=axDiffVx)
        axDiffVy.tick_params('y',labelleft = False)
    else:
        axDiffVy = plt.subplot(335)
        axDiffVy.set_ylabel("Velocity (km/s)")
    axDiffVy.tick_params('both',length = 12)
    axDiffVy.plot(orbit0.t,diff_Vy)
    axDiffVy.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVy.set_title("Difference in $v_y$")
    axDiffVy.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Vy,diff_Vz, 0.1):
        axDiffVz = plt.subplot(336,sharex = axDiffZ, sharey=axDiffVy)
        axDiffVz.tick_params('y',labelleft = False)
    else:
        axDiffVz = plt.subplot(336)
        axDiffVz.set_ylabel("Velocity (km/s)")
    axDiffVz.plot(orbit0.t,diff_Vz)
    axDiffVz.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVz.set_title("Difference in $v_z$")
    axDiffVz.tick_params('x',labelbottom = False)
    axDiffVz.tick_params('both',length = 12)
    
    magPosDifference = np.sqrt(diff_X**2+diff_Y**2+diff_Z**2)
    
    axDiffMagPos = plt.subplot(337,sharex = axDiffVx)
    axDiffMagPos.plot(orbit0.t,magPosDifference)
    axDiffMagPos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagPos.set_title("Magnitude of difference in position")
    axDiffMagPos.set_ylabel(r"$|\vec x|$ (kpc)")
    axDiffMagPos.set_xlabel("t (Myr)")
    axDiffMagPos.tick_params('both',length = 12)
    
    magVelDifference = np.sqrt(diff_Vx**2+diff_Vy**2+diff_Vz**2)
    
    axDiffMagVel = plt.subplot(339,sharex = axDiffVz)
    axDiffMagVel.plot(orbit0.t,magVelDifference)
    axDiffMagVel.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagVel.set_title("Magnitude of difference in velocity")
    axDiffMagVel.set_ylabel(r"$|\vec v|$ (km/s)")
    axDiffMagVel.set_xlabel("t (Myr)")
    axDiffMagVel.tick_params('both',length = 12)
    
    mean_Vx = 0.5*(orbit1.v_x+orbit0.v_x).to(u.km/u.s)
    mean_Vy = 0.5*(orbit1.v_y+orbit0.v_y).to(u.km/u.s)
    mean_Vz = 0.5*(orbit1.v_z+orbit0.v_z).to(u.km/u.s)
    
    magVelMean = np.sqrt(mean_Vx**2+mean_Vy**2+mean_Vz**2)
    
    dirCos = []
    for (i, pos) in enumerate(diff_X):
        PosDotVel = 0.5*diff_X[i]*mean_Vx[i]+0.5*diff_Y[i]*mean_Vy[i]+0.5*diff_Z[i]*mean_Vz[i]
        dirCos.append(PosDotVel / (0.5*magPosDifference[i] * magVelMean[i]))
    
    axDirCos = plt.subplot(338,sharex = axDiffVy)
    axDirCos.plot(orbit0.t,dirCos)
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$",fontsize = 42)
    axDirCos.set_xlabel("t (Myr)")
    axDirCos.tick_params('both',length = 12)
    
    diffPath = dirPathPair+"/VariableDifference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
    if(saveGraphs): 
       plt.savefig(diffPath)
       
    plt.close(fig1)
    plt.close(fig2)
    
    
    
