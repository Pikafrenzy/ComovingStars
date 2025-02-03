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

#initial variables
mw = gp.MilkyWayPotential()
T = 1000*u.Myr
def trunc3dp(number):
    return np.trunc(number*1e3)/1e3
    
def makeLabel(pos, vel):
    label = "Position = "
    label += "["+str(trunc3dp(pos[0].to_value()))+", "+str(trunc3dp(pos[1].to_value()))+", "+str(trunc3dp(pos[2].to_value()))+"] "+pos[0].unit.to_string()
    label += ", Velocity = "
    label += "["+str(trunc3dp(vel[0].to_value()))+", "+str(trunc3dp(vel[1].to_value()))+", "+str(trunc3dp(vel[2].to_value()))+"] "+vel[0].unit.to_string()
    return label

def pairGraph(ID,star0, star1, saveGraphs, dirPath):
    w0 = gd.PhaseSpacePosition(pos = star0.get_Pos(),vel = star0.get_Vel())
    w1 = gd.PhaseSpacePosition(pos = star1.get_Pos(),vel = star1.get_Vel())
    orbit0 = mw.integrate_orbit(w0, dt = 1*u.Myr, t1=0, t2 = T)
    orbit1 = mw.integrate_orbit(w1, dt = 1*u.Myr, t1=0, t2 = T)
    
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(16, 12),layout = 'constrained')
    gs = gridspec.GridSpec(4, 12, figure=fig, height_ratios=[2, 1, 1, 1],wspace=0.2,hspace=0.2)  
    axY = fig.add_subplot(gs[0,0:3])
    axZ = fig.add_subplot(gs[0,3:6],sharey = axY)
    axVelY = fig.add_subplot(gs[0,6:9])
    axVelZ = fig.add_subplot(gs[0,9:12])
    
    axY.plot(orbit0.pos.x,orbit0.pos.y,label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axY.plot(orbit1.pos.x,orbit1.pos.y,label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axY.set_box_aspect(1)
    axY.axis('equal')
    axY.set_title("Y against X")
    axY.set_ylabel("y (kpc), z (kpc)")
    axY.set_xlabel("x (kpc)")
    
    axZ.plot(orbit0.pos.x,orbit0.pos.z,label=makeLabel(star0.get_Pos(), star0.get_Vel()),linewidth = 0.3, color = 'r')
    axZ.plot(orbit1.pos.x,orbit1.pos.z,label=makeLabel(star1.get_Pos(), star1.get_Vel()),linewidth = 0.3, color = 'b')
    axZ.set_box_aspect(1)
    axZ.axis('equal')
    axZ.legend(fontsize=10, bbox_to_anchor=(1,1.3))
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
    
    diff_X = orbit1.pos.x-orbit0.pos.x
    diff_Y = orbit1.pos.y-orbit0.pos.y
    diff_Z = orbit1.pos.z-orbit0.pos.z
    diff_Vx = (orbit1.v_x-orbit0.v_x).to(u.km/u.s)
    diff_Vy = (orbit1.v_y-orbit0.v_y).to(u.km/u.s)
    diff_Vz = (orbit1.v_z-orbit0.v_z).to(u.km/u.s)
        
    #position plots
    axDiffX = fig.add_subplot(gs[1,0:4])
    axDiffX.plot(orbit0.t,diff_X)
    axDiffX.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffX.set_title("Difference in x")
    axDiffX.set_ylabel("Position (kpc)")
    axDiffX.tick_params('x',labelbottom = False)
    
    def checkSharing(x1, x2, margin):
        diff_min = x2.min()-x1.min()
        diff_max = x2.max()-x1.max()
        share = diff_min.to_value() <= margin and diff_max.to_value() <= margin
        return share
    
    if checkSharing(diff_X,diff_Y,0.01):
        axDiffY = fig.add_subplot(gs[1,4:8],sharey=axDiffX)
        axDiffY.tick_params('y',labelleft = False)
    else:
        axDiffY = fig.add_subplot(gs[1,4:8])
        axDiffY.set_ylabel("Position (kpc)")
    axDiffY.plot(orbit0.t,diff_Y)
    axDiffY.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffY.set_title("Difference in y")
    axDiffY.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Y,diff_Z, 0.01):
        axDiffZ = fig.add_subplot(gs[1,8:12],sharey=axDiffY)
        axDiffZ.tick_params('y',labelleft = False)
    else:
        axDiffZ = fig.add_subplot(gs[1,8:12])
        axDiffZ.set_ylabel("Position (kpc)")
    axDiffZ.plot(orbit0.t,diff_Z)
    axDiffZ.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffZ.set_title("Difference in z")
    axDiffZ.tick_params('x',labelbottom = False)
    
    #velocity plots
    axDiffVx = fig.add_subplot(gs[2,0:4],sharex=axDiffX)
    axDiffVx.plot(orbit0.t,diff_Vx)
    axDiffVx.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVx.set_title("Difference in $v_x$")
    axDiffVx.set_ylabel("Velocity (km/s)")
    axDiffVx.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Vx,diff_Vy, 0.1):
        axDiffVy = fig.add_subplot(gs[2,4:8],sharex = axDiffY, sharey=axDiffVx)
        axDiffVy.tick_params('y',labelleft = False)
    else:
        axDiffVy = fig.add_subplot(gs[2,4:8])
        axDiffVy.set_ylabel("Velocity (km/s)")
    axDiffVy.plot(orbit0.t,diff_Vy)
    axDiffVy.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVy.set_title("Difference in $v_y$")
    axDiffVy.tick_params('x',labelbottom = False)
    
    if checkSharing(diff_Vy,diff_Vz, 0.1):
        axDiffVz = fig.add_subplot(gs[2,8:12],sharex = axDiffZ, sharey=axDiffVy)
        axDiffVz.tick_params('y',labelleft = False)
    else:
        axDiffVz = fig.add_subplot(gs[2,8:12])
        axDiffVz.set_ylabel("Velocity (km/s)")
    axDiffVz.plot(orbit0.t,diff_Vz)
    axDiffVz.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffVz.set_title("Difference in $v_z$")
    axDiffVz.tick_params('x',labelbottom = False)
    
    magPosDifference = np.sqrt(diff_X**2+diff_Y**2+diff_Z**2)
    
    axDiffMagPos = fig.add_subplot(gs[3,0:4],sharex = axDiffVx)
    axDiffMagPos.plot(orbit0.t,magPosDifference)
    axDiffMagPos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagPos.set_title("Magnitude of difference in position")
    axDiffMagPos.set_ylabel(r"$|\vec x|$ (kpc)")
    axDiffMagPos.set_xlabel("t (Myr)")
    
    magVelDifference = np.sqrt(diff_Vx**2+diff_Vy**2+diff_Vz**2)
    
    axDiffMagVel = fig.add_subplot(gs[3,8:12],sharex = axDiffVz)
    axDiffMagVel.plot(orbit0.t,magVelDifference)
    axDiffMagVel.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDiffMagVel.set_title("Magnitude of difference in velocity")
    axDiffMagVel.set_ylabel(r"$|\vec v|$ (km/s)")
    axDiffMagVel.set_xlabel("t (Myr)")
    
    mean_Vx = 0.5*(orbit1.v_x+orbit0.v_x).to(u.km/u.s)
    mean_Vy = 0.5*(orbit1.v_y+orbit0.v_y).to(u.km/u.s)
    mean_Vz = 0.5*(orbit1.v_z+orbit0.v_z).to(u.km/u.s)
    
    magVelMean = np.sqrt(mean_Vx**2+mean_Vy**2+mean_Vz**2)
    
    dirCos = []
    for (i, pos) in enumerate(diff_X):
        PosDotVel = 0.5*diff_X[i]*mean_Vx[i]+0.5*diff_Y[i]*mean_Vy[i]+0.5*diff_Z[i]*mean_Vz[i]
        dirCos.append(PosDotVel / (0.5*magPosDifference[i] * magVelMean[i]))
    
    axDirCos = fig.add_subplot(gs[3,4:8],sharex = axDiffVy)
    axDirCos.plot(orbit0.t,dirCos)    
    axDirCos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
    axDirCos.set_title("Direction Cosine")
    axDirCos.set_ylabel(r"$\frac{x \cdot \overline{v}}{|x||\overline{v}|}$")
    axDirCos.set_xlabel("t (Myr)")
    
    pathPair = dirPath+"/Star_Pair_"+str(ID)+".png"
    if(saveGraphs): 
        plt.savefig(pathPair)   
    else:
        plt.show()
    plt.close(fig)
    
    
    
