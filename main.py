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

mw = gp.MilkyWayPotential()
T = 500*u.Myr
pos0 = [8.0,0.0,2.0]*u.kpc
posDiff = [0.0,0.0,0.001]*u.kpc
pos1= pos0 + posDiff
vel0 = [20.0,200.0,30.0]*u.km/u.s
velDiff = [0.0,0.0,0.0]*u.km/u.s
vel1 = vel0 + velDiff

saveGraphs = False

w0 = gd.PhaseSpacePosition(pos=pos0, vel=vel0)
w1 = gd.PhaseSpacePosition(pos=pos1, vel=vel1)
orbit0 = mw.integrate_orbit(w0, dt = 1*u.Myr, t1=0, t2=T)
orbit1 = mw.integrate_orbit(w1, dt = 1*u.Myr, t1=0, t2=T)

plt.rcParams.update({'font.size': 10})

fig1, ax1 = plt.subplots(1,2, layout="constrained")
ax1[0].plot(orbit0.pos.x,orbit0.pos.y,label="Position = [8.0,0.0,2.0]",linewidth = 0.3, color = 'r')
ax1[0].plot(orbit1.pos.x,orbit1.pos.y,label="Position = [8.0,0.0,2.001]",linewidth = 0.3, color = 'b')
ax1[0].legend(fontsize=10, bbox_to_anchor=(0.7,-0.2))
ax1[0].set_aspect('equal','box')
ax1[0].set_title("Y against X")

ax1[1].plot(orbit0.pos.x,orbit0.pos.z,label="Position = [8.0,0.0,2.0]",linewidth = 0.3, color = 'r')
ax1[1].plot(orbit1.pos.x,orbit1.pos.z,label="Position = [8.0,0.0,2.001]",linewidth = 0.3, color = 'b')
ax1[1].axis('equal')
ax1[1].set(ylim=(-4,4))
ax1[1].set_title("Z against X")

dirPath = "../Plots/"+datetime.now().strftime("%Y%m%d_%H%M%S")
if(saveGraphs): 
   os.mkdir(dirPath)
posPath = dirPath+"/Position_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
    plt.savefig(posPath)

layout = [['DiffX','DiffY','DiffZ'],
          ['DiffVX','DiffVY','DiffVZ'],
          ['DiffPos','.','DiffVel']]

plt.rcParams.update({'font.size': 28})

fig2, ax2 = plt.subplot_mosaic(layout, figsize=(30, 10),layout="constrained")

ax2['DiffX'].plot(orbit0.t,orbit1.pos.x-orbit0.pos.x)
ax2['DiffX'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffX'].set_title("Difference in x")
ax2['DiffX'].set_ylabel("x (kpc)")
ax2['DiffX'].set_xlabel("t (Myr)")

ax2['DiffY'].plot(orbit0.t,orbit1.pos.y-orbit0.pos.y)
ax2['DiffY'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffY'].set_title("Difference in y")
ax2['DiffY'].set_ylabel("y (kpc)")
ax2['DiffY'].set_xlabel("t (Myr)")

ax2['DiffZ'].plot(orbit0.t,orbit1.pos.z-orbit0.pos.z)
ax2['DiffZ'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffZ'].set_title("Difference in z")
ax2['DiffZ'].set_ylabel("z (kpc)")
ax2['DiffZ'].set_xlabel("t (Myr)")

ax2['DiffVX'].plot(orbit0.t,(orbit1.v_x-orbit0.v_x).to(u.km/u.s))
ax2['DiffVX'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffVX'].set_title("Difference in $v_x$")
ax2['DiffVX'].set_ylabel("$v_x$ (km/s)")
ax2['DiffVX'].set_xlabel("t (Myr)")

ax2['DiffVY'].plot(orbit0.t,(orbit1.v_y-orbit0.v_y).to(u.km/u.s))
ax2['DiffVY'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffVY'].set_title("Difference in $v_y$")
ax2['DiffVY'].set_ylabel("$v_y$ (km/s)")
ax2['DiffVY'].set_xlabel("t (Myr)")

ax2['DiffVZ'].plot(orbit0.t,(orbit1.v_z-orbit0.v_z).to(u.km/u.s))
ax2['DiffVZ'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffVZ'].set_title("Difference in $v_z$")
ax2['DiffVZ'].set_ylabel("$v_z$ (km/s)")
ax2['DiffVZ'].set_xlabel("t (Myr)")

magPos_orbit0 = []
for (i,pos) in enumerate(orbit0.pos.x):
    magPos_i = np.sqrt(orbit0.pos.x[i]**2+orbit0.pos.y[i]**2+orbit0.pos.z[i]**2)
    magPos_orbit0.append(magPos_i)

magPos_orbit1 = []
for (j,pos) in enumerate(orbit1.pos.x):
    magPos_j = np.sqrt(orbit1.pos.x[j]**2+orbit1.pos.y[j]**2+orbit1.pos.z[j]**2)
    magPos_orbit1.append(magPos_j)

magPosDifference = []
for (k,pos) in enumerate(magPos_orbit0):
    magPosDifference.append(magPos_orbit1[k].value-magPos_orbit0[k].value)

ax2['DiffPos'].plot(orbit0.t,magPosDifference)
ax2['DiffPos'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffPos'].set_title("Difference in magnitude of position")
ax2['DiffPos'].set_ylabel(r"$|\vec x|$ (kpc)")
ax2['DiffPos'].set_xlabel("t (Myr)")

magVel_orbit0 = []
for (i,vel) in enumerate(orbit0.v_x):
    magVel_i = np.sqrt(orbit0.v_x[i]**2+orbit0.v_y[i]**2+orbit0.v_z[i]**2)
    magVel_orbit0.append(magVel_i)

magVel_orbit1 = []
for (j,vel) in enumerate(orbit1.v_x):
    magVel_j = np.sqrt(orbit1.v_x[j]**2+orbit1.v_y[j]**2+orbit1.v_z[j]**2)
    magVel_orbit1.append(magVel_j)

magVelDifference = []
for (k,pos) in enumerate(magVel_orbit0):
    magVelDifference.append(magVel_orbit1[k].to(u.km/u.s).value-magVel_orbit0[k].to(u.km/u.s).value)

ax2['DiffVel'].plot(orbit0.t,magVelDifference)
ax2['DiffVel'].plot(orbit0.t, orbit0.t*0, color = (0.0,0.0,0.0,0.5))
ax2['DiffVel'].set_title("Difference in magnitude of velocity")
ax2['DiffVel'].set_ylabel(r"$|\vec v|$ (km/s)")
ax2['DiffVel'].set_xlabel("t (Myr)")

diffPath = dirPath+"/VariableDifference_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".png"
if(saveGraphs): 
   plt.savefig(diffPath)




