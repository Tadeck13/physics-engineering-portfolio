"""
Script that generates circles offest from the centre to avoid the region of extreme sensitivity.
The cricle motion have defineable centres, radii, periods etc etc.  

Author: Rayhaan Perin & Tadeck Jones
Date: 03-09-2024
Packages:
Optional Packages:
Changelist:  
            -- 03-09-2024 -- Creation of script.
"""

import numpy as np
from tqdm import tqdm 
import csv

def circleMotion(deltaT: float, radius: float, tangentialVeloicty: float, offset: float) -> np.ndarray:
    """
    Defintiion of cricular motion offset from centre.  Period 
    is calcualted from radius and tangential veloicty.  

    Input:
            deltaT (float) - The change in time/time step
            radius (float) - The radius of the cricle
            tangentialVelocity (float) - The tangential velocity of the circle
            offset (float) - The offset in the x, y and z dimensions
    Output:
            path (np.ndarray) - The (t (s), x (mm), y (mm), z (mm)) path where the parameters are columns.  
    
    """
    # Conversion from radius and tangential velocity to period
    # This is from w = 2pi/T and v = w x r
    period = (2*np.pi*radius)/tangentialVeloicty

    # time , simulate for twice the period
    t = np.arange(0, 2*period + deltaT, deltaT)

    # positions 
    x = radius*np.cos(((2*np.pi)/period)*t) + offset
    y = radius*np.sin(((2*np.pi)/period)*t) + offset
    z = np.repeat(offset, len(x))

    # path
    path = np.array([t, x, y, z], dtype = np.float32).T

    return path

def writeData(data: np.ndarray, OUTPATH: str, OUTNAME: str) -> None:
    """
    Writing of the *.placements file is performed by this function.  

    Input:
            data (np.ndarray) - The input data to write.  The columns are (t, x, y, z) in s and mm respectively.  
            OUTPATH (str) - The output path/folder.  
            OUTNAME (str) - The output file name.  
    """
    # define zeros and ones for placement columns that are irrelevant 
    zero = np.array([0 for _ in range(len(data[:, 0]))], dtype = np.int64)
    one = np.array([1 for _ in range(len(data[:, 0]))], dtype = np.int64)

    # The array to write
    writeData = np.array([data[:, 0], zero, zero, one, zero, data[:, 1], data[:, 2], data[:, 3]]).T # The addition - min is to start the motion at 0.0 s within FP precision

    # Add header and write data
    with open("{}{}".format(OUTPATH, OUTNAME), 'w') as f:
        f.write("Time s\n")
        f.write("Rotation deg\n")
        f.write("Translation mm\n")
        csv.writer(f, delimiter=' ').writerows(writeData)

def main():
    OUTPATH = "/home/rayhaan/REPO_HR++/GATE_HR/Circles/Placements/"
    radii = [5, 4, 3, 2, 1, 0.5, 0.1]
    deltaT = 0.00001 
    tangnetialVelocity = 10.0

    for i in tqdm(range(len(radii))):
        # Create path
        path = circleMotion(deltaT = deltaT, radius = radii[i], tangentialVeloicty = tangnetialVelocity, offset = 50)

        # Define output file name 
        OUTNAME = "cirlce_vT_{}mps_t_2period_deltaT_{}s_radius_{}mm.placements".format(tangnetialVelocity, deltaT, radii[i])

        # Write the data
        writeData(data = path, OUTPATH = OUTPATH, OUTNAME = OUTNAME)

if __name__ == "__main__":
    main()