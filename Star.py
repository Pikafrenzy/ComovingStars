#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 15:08:20 2024

@author: ambroselo
"""

class Star:
    def __init__(self, posx,posy,posz,vx,vy,vz):
        self.Position = [posx,posy,posz]
        self.Velocity = [vx,vy,vz]
        
    def __str__(self):
        output = "Position: ["
        for i in self.Position[:-1]:
            output += str(i) + ", "
        output += str(self.Position[-1]) +"] \n Velocity: ["
        for j in self.Velocity[:-1]:
            output += str(j) + ", "
        output += str(self.Velocity[-1]) +"]"
        return output