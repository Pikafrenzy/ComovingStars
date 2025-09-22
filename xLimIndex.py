#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 08:59:56 2025

@author: ambroselo
"""

import astropy.units as u

def xLimIndex(startTime, limit, dt):
    steps = (limit - startTime)/dt
    return steps