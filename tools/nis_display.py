# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def nis_ana(file_name):
    df = pd.read_table(file_name)
    
    #print(df.iloc[[2]])
    
    sensor_types = ['radar', 'lidar']
    thresholds = [7.815, 5.991]
    for sensor,th in zip(sensor_types, thresholds):
        radar_rows=df.loc[df['sensor_type']==sensor]
        nis=radar_rows['NIS']
    
        stat = np.mean(nis>th)
    
        print(stat)
        
        plt.figure(figsize=(12,9))
        plt.plot(nis)
        plt.axhline(y=th)
        plt.title('{} NIS={}'.format(sensor, stat))
        plt.ylim((0,20))
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='NIS analysis')
    parser.add_argument('file',type=str,help='output file from UnscentedKF',
                        default='result.txt')
    
    args = parser.parse_args()
    
    file = args.file
    
    nis_ana(file)
    