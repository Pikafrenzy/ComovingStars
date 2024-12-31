#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 15:08:20 2024

@author: ambroselo
"""

from collections import namedtuple

class Star:
    def __init__(self, posx,posy,posz,v_x,v_y,v_z):
        self.Position = namedtuple('Position',['x','y','z'])
        
        self.Position.x = posx
        self.Position.y = posy
        self.Position.z = posz
        
        self.Velocity = namedtuple('Position',['v_x','v_y','v_z'])
        
        self.Velocity.v_x = v_x
        self.Velocity.v_y = v_y
        self.Velocity.v_z = v_z
        
    def __str__(self):
        output = "Position: ["
        for i in self.Position[:-1]:
            output += str(i) + ", "
        output += str(self.Position[-1]) +"] \n Velocity: ["
        for j in self.Velocity[:-1]:
            output += str(j) + ", "
        output += str(self.Velocity[-1]) +"]"
        return output
    
    def get_x(self):
        return self.Position.x
    def get_y(self):
        return self.Position.y
    def get_z(self):
        return self.Position.z
    def get_Vx(self):
        return self.Velocity.v_x
    def get_Vy(self):
        return self.Velocity.v_y
    def get_Vz(self):
        return self.Velocity.v_z