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

# Gala
import gala.dynamics as gd
import gala.potential as gp
import gala.units as gu

#initial variables
mw = gp.MilkyWayPotential()
T = 500*u.Myr
pos0 = [8.0,0.0,2.0]*u.kpc
posDiff = [0.0,0.0,0.001]*u.kpc
pos1 = pos0 + posDiff
vel0 = [20.0,200.0,30.0]*u.km/u.s
velDiff = [0.0,0.0,0.0]*u.km/u.s
vel1 = vel0 + velDiff

#toggle for whether graphs are saved automatically to computer
saveGraphs = False

w0 = gd.PhaseSpacePosition(pos=pos0, vel=vel0)
w1 = gd.PhaseSpacePosition(pos=pos1, vel=vel1)
orbit0 = mw.integrate_orbit(w0, dt = 1*u.Myr, t1=0, t2=T)
orbit1 = mw.integrate_orbit(w1, dt = 1*u.Myr, t1=0, t2=T)

plt.rcParams.update({'font.size': 10})

def trunc3dp(number):
    return np.trunc(number*1e3)/1e3
    
def makeLabel(vec3):
    label = "["+str(trunc3dp(vec3[0].to_value()))+", "+str(trunc3dp(vec3[1].to_value()))+", "+str(trunc3dp(vec3[2].to_value()))+"] "+vec3[0].unit.to_string()
    return label

fig1, ax1 = plt.subplots(1,2, layout="constrained")
ax1[0].plot(orbit0.pos.x,orbit0.pos.y,label="Position = "+makeLabel(pos0),linewidth = 0.3, color = 'r')
ax1[0].plot(orbit1.pos.x,orbit1.pos.y,label="Position = "+makeLabel(pos1),linewidth = 0.3, color = 'b')
ax1[0].legend(fontsize=10, bbox_to_anchor=(0.9,-0.2))
ax1[0].set_aspect('equal','box')
ax1[0].set_title("Y against X")

ax1[1].plot(orbit0.pos.x,orbit0.pos.z,label="Position = "+makeLabel(pos0),linewidth = 0.3, color = 'r')
ax1[1].plot(orbit1.pos.x,orbit1.pos.z,label="Position = "+makeLabel(pos1),linewidth = 0.3, color = 'b')
ax1[1].axis('equal')
ax1[1].set(ylim=(-4,4))
ax1[1].set_title("Z against X")

dirPath = "../Plots/"+datetime.now().strftime("%Y%m%d_%H%M%S")
if(saveGraphs): 
   os.mkdir(dirPath)
posPath = dirPath+"/Position_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
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
axDiffVy.set_xlabel("t (Myr)")

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


magPos_orbitx = []
for (i,posx) in enumerate(orbit0.pos.x):
    magPos_i = orbit1.pos.x[i]-orbit0.pos.x[i]
    magPos_orbitx.append(magPos_i)

magPos_orbity = []
for (j,posy) in enumerate(orbit0.pos.y):
    magPos_j = orbit1.pos.y[j]-orbit0.pos.y[j]
    magPos_orbity.append(magPos_j)

magPos_orbitz = []
for (k,posz) in enumerate(orbit0.pos.z):
    magPos_k = orbit1.pos.z[k]-orbit0.pos.z[k]
    magPos_orbitz.append(magPos_k)
    
magPosDifference = []
for (l,pos) in enumerate(magPos_orbitx):
    magPosDifference.append(np.sqrt(magPos_orbitx[l].value**2+magPos_orbity[l].value**2+magPos_orbitz[l].value**2))

axDiffMagPos = plt.subplot(337,sharex = axDiffVx)
axDiffMagPos.plot(orbit0.t,magPosDifference)
axDiffMagPos.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
axDiffMagPos.set_title("Magnitude of the difference in position")
axDiffMagPos.set_ylabel(r"$|\vec x|$ (kpc)")
axDiffMagPos.set_xlabel("t (Myr)")
axDiffMagPos.tick_params('both',length = 12)

magVel_orbitx = []
for (i,velx) in enumerate(orbit0.v_x):
    magVel_i = (orbit1.v_x[i]-orbit0.v_x[i]).to(u.km/u.s)
    magVel_orbitx.append(magVel_i)

magVel_orbity = []
for (j,vely) in enumerate(orbit0.v_y):
    magVel_j = (orbit1.v_y[j]-orbit0.v_y[j]).to(u.km/u.s)
    magVel_orbity.append(magVel_j)
    
magVel_orbitz = []
for (k,velz) in enumerate(orbit0.v_z):
    magVel_k = (orbit1.v_z[k]-orbit0.v_z[k]).to(u.km/u.s)
    magVel_orbitz.append(magVel_k)

magVelDifference = []
for (l,vel) in enumerate(magVel_orbitx):
    magVelDifference.append(np.sqrt(magVel_orbitx[l].value**2+magVel_orbity[l].value**2+magVel_orbitz[l].value**2))
    
axDiffMagVel = plt.subplot(339,sharex = axDiffVz)
axDiffMagVel.plot(orbit0.t,magVelDifference)
axDiffMagVel.plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
axDiffMagVel.set_title("DMagnitude of the difference in velocity")
axDiffMagVel.set_ylabel(r"$|\vec v|$ (km/s)")
axDiffMagVel.set_xlabel("t (Myr)")
axDiffMagVel.tick_params('both',length = 12)

diffPath = dirPath+"/VariableDifference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
   plt.savefig(diffPath)



