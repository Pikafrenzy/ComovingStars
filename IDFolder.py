#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 22:17:51 2025

@author: ambroselo
"""
import os

def IDdirCheck(ID, dirPath, saveGraphs):
    IDdirPath = dirPath/str(ID)
    if not IDdirPath.is_dir(): 
        if saveGraphs:
            os.mkdir(IDdirPath)
    return IDdirPath