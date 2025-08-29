#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 15:21:43 2025

@author: ambroselo
"""
import numpy as np
from datetime import datetime
import time
import initConditionGenerator as ICG

# notes the start time of running the program to determine how long it ran for
startTime = time.time()

# toggle for whether graphs are saved automatically to computer
saveGraphs = False
# toggle for whether stars are integrated
integrateOrbits = False

# saves start time for directories
dirTime = datetime.now().strftime("%Y%m%d_%H%M%S")

# seeded RNGs to ensure consistent data when running multiple times
# two RNGs ensures adding more stars doesn't change the original stars
centrePosRNG = np.random.default_rng(137) 
differenceRNG = np.random.default_rng(1836)

starPairCount = 10000
pairGraphLimit = 12

# fiducial point initial distribution graphs
starCentres, starCoords = ICG.generateStarCentres(starPairCount, centrePosRNG)
ICG.starCentresGraph(starCoords, saveGraphs,dirTime)

# pair graphs
starPairs = ICG.starPairsCreate(differenceRNG, starCentres)
ICG.starPairGraphs(starPairs, saveGraphs,integrateOrbits, pairGraphLimit, dirTime)

# graphing invariants of initial positions
ICG.invariantGraphs(*ICG.invariantArrays(starPairs),saveGraphs, dirTime)

# outputs the time it took the program to run
endTime = time.time()
print("This took " + str(endTime-startTime) + " seconds")