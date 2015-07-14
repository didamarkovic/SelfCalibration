#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes / David Hogg
# August 2012

# Modified by Dida Markovic
# 2015

# This script is used to generate the four survey strategies investigated in
# the self-calibration paper.

from __future__ import division, print_function

import numpy as np
import parameters


def uniform_center_list(sky_limits, FoV, nx, ny):
    assert nx * FoV[0] > sky_limits[1] - sky_limits[0]
    assert ny * FoV[1] > sky_limits[3] - sky_limits[2]
    
    x_step_size = (sky_limits[1] - sky_limits[0] - FoV[0]) / (nx - 1.)
    y_step_size = (sky_limits[3] - sky_limits[2] - FoV[1]) / (ny - 1.)

    x_center_list = sky_limits[0] + 0.5 * FoV[0] + x_step_size * np.arange(nx)
    y_center_list = sky_limits[2] + 0.5 * FoV[1] + y_step_size * np.arange(ny)
    return x_center_list, y_center_list


def generate_uniform_survey(sky_limits, FoV, Ncovering,
                                                rotate=False, offset=False):

    theta = 0.
    nx = 12
    ny = 12

    x_center_list, y_center_list = uniform_center_list(sky_limits, FoV, nx, ny)
    if offset:
        x_center_list_1, y_center_list_1 =\
                    uniform_center_list(sky_limits, FoV, nx + 1, ny - 1)
        x_center_list_2, y_center_list_2 =\
                    uniform_center_list(sky_limits, FoV, nx - 1, ny + 1)

    x = np.zeros((Ncovering * nx * ny, 4))
    indx = 0
    for covering in range(Ncovering):
        xcl = x_center_list
        ycl = y_center_list
        tnx = nx
        tny = ny
        if offset:
            if covering % 3 == 1:
                xcl = x_center_list_1
                ycl = y_center_list_1
                tnx = nx + 1
                tny = ny - 1
            if covering % 3 == 2:
                xcl = x_center_list_2
                ycl = y_center_list_2
                tnx = nx - 1
                tny = ny + 1
        for yy in range(tny):
            x[indx:indx + tnx, 0] = indx + np.arange(tnx)
            x[indx:indx + tnx, 1] = xcl
            x[indx:indx + tnx, 2] = ycl[yy]
            x[indx:indx + tnx, 3] = theta
            indx += tnx
        if rotate:
            theta += 360. / Ncovering

    return x[:indx, :]


def generate_random_survey(srvy, FoV):
    N = len(srvy[:, 0])    
    # Add random offsets
    srvy[:, 1] += np.random.uniform(-0.5 * FoV[0], 0.5 * FoV[0], size=N)
    srvy[:, 2] += np.random.uniform(-0.5 * FoV[1], 0.5 * FoV[1], size=N)
    srvy[:, 3] += np.random.uniform(0., 360, size=N)
    return srvy

def generate_euclid_survey(sky_limits, FoV, pattern, rand=False):

    nx = ny = 12
    Ncovering = 4
    rotate = theta = 0.0

    x_center_list, y_center_list = uniform_center_list(sky_limits, FoV, nx, ny)

    x = np.zeros((Ncovering * nx * ny, 4))
    indx = 0
    displace = [0.0, 0.0]
    for covering in range(Ncovering):
        xcl = x_center_list
        ycl = y_center_list
        tnx = nx
        tny = ny

        if sum(sum(pattern)) > 0.0:
            idith = covering % 4
            if idith > 0:
                displace += pattern[idith-1]
                xcl = xcl + displace[0]
                ycl = ycl + displace[1]

        for yy in range(tny):
            x[indx:indx + tnx, 0] = indx + np.arange(tnx)
            x[indx:indx + tnx, 1] = xcl
            x[indx:indx + tnx, 2] = ycl[yy]
            x[indx:indx + tnx, 3] = theta
            indx += tnx
        if rotate:
            theta += 360. / Ncovering

    return x[:indx, :]

if __name__ == "__main__":

    dic = eval(open('parameters.py').read())
    sky_limits = [-4., 4., -4., 4.]
    number_passes = 4
    FoV = dic['FoV']

    xC = generate_uniform_survey(sky_limits, FoV, number_passes, offset=True)
    np.savetxt('C.txt', xC)

    xD = generate_random_survey(xC, FoV)
    np.savetxt('D.txt', xD)

    # Euclid
    randomise=0.0 # How much to randomise the dither steps (TODO)
    Euclid_FoV = [4*612. + 3*50., 4*612. + 3*100.] # arcsec

    J = np.array([[50.0, 100.0], [0.0, 100.0], [0.0, 100.0]])
    xJ = generate_euclid_survey(sky_limits, FoV, pattern=J/Euclid_FoV*FoV, rand=randomise)
    np.savetxt('J.txt', xJ)

    S = np.array([[50.0, 100.0], [0.0, 100.0], [50.0, 100.0]])
    xS = generate_euclid_survey(sky_limits, FoV, pattern=S/Euclid_FoV*FoV, rand=randomise)
    np.savetxt('S.txt', xS)
