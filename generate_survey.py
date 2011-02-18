#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import random
import tparameters as p
import matplotlib.pylab as plt

# XX Dangerous hack!!
pdic = p.dic
FoV = pdic['FoV']
sky_limits = pdic['sky_limits']

def fp2sky(x, y, pointing, orientation):
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
  beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
  return alpha, beta

def plot_survey(survey, survey_name):
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  x = np.array([x_min, x_min, x_max, x_max, x_min])
  y = np.array([y_min, y_max, y_max, y_min, y_min])
  for image in survey:
    alpha, beta = fp2sky(x,y,image[1:3], image[3])
    plt.plot(alpha,beta,'k-',alpha=0.25)
  plt.xlabel(r"$\alpha$", fontsize=25)
  plt.ylabel(r"$\beta$", fontsize=25)
  plt.title(r"Strategy %s: %i Pointings" % (survey_name, len(survey[:,0])))
  plt.axis('equal')
  plt.xlim(sky_limits[0]-FoV[0], sky_limits[1]+FoV[0])
  plt.ylim(sky_limits[2]-FoV[1], sky_limits[3]+FoV[1])
  return

def uniform_center_list(nx,ny):
  x_step_size = (sky_limits[1]-sky_limits[0] - FoV[0])/(nx-1.)
  y_step_size = (sky_limits[3]-sky_limits[2] - FoV[1])/(ny-1.)
  
  x_center_list = sky_limits[0] + 0.5*FoV[0] + x_step_size*np.arange(nx)
  y_center_list = sky_limits[2] + 0.5*FoV[1] + y_step_size*np.arange(ny)
  print 'uniform_center_list: overlap: %.2f %% %.2f %% ' % ((100*((FoV[0]-x_step_size)/FoV[0])), (100*((FoV[1]-y_step_size)/FoV[1])))
  return (x_center_list, y_center_list)

def generate_uniform_survey(Ncovering, rotate = False, offset = False):
  theta = 0.
  # XX why 5 %?
  nx = np.ceil((sky_limits[1]-sky_limits[0])/(0.975*FoV[0])).astype(int)
  ny = np.ceil((sky_limits[3]-sky_limits[2])/(0.975*FoV[1])).astype(int)
  (x_center_list, y_center_list) = uniform_center_list(nx,ny)
  if offset:
    (x_center_list_1, y_center_list_1) = uniform_center_list(nx+1,ny-1)
    (x_center_list_2, y_center_list_2) = uniform_center_list(nx-1,ny+1)
  
  x = np.zeros((Ncovering*nx*ny,4))
  ii = 0
  for covering in range(Ncovering):
    xcl = x_center_list
    ycl = y_center_list
    tnx = nx
    tny = ny
    if offset:
      if covering%3 == 1:
        xcl = x_center_list_1
        ycl = y_center_list_1
        tnx = nx+1
        tny = ny -1
      if covering%3 == 2:
        xcl = x_center_list_2
        ycl = y_center_list_2
        tnx = nx-1
        tny = ny+1
    for yy in range(tny):
      x[ii:ii+tnx,0] = ii+np.arange(tnx)
      x[ii:ii+tnx,1] = xcl
      x[ii:ii+tnx,2] = ycl[yy]
      x[ii:ii+tnx,3] = theta
      ii += tnx
    if rotate:
      theta += 360./Ncovering 
  return x[:ii,:]

def generate_random_survey(N):
  nx = np.ceil((sky_limits[1]-sky_limits[0])/FoV[0]).astype(int)
  ny = np.ceil((sky_limits[3]-sky_limits[2])/FoV[1]).astype(int)
  x_box_size = (sky_limits[1] - sky_limits[0])/nx
  y_box_size = (sky_limits[3] - sky_limits[2])/ny
  ntheta = np.floor(N/(nx*ny)).astype(int)
  theta_box_size = 360./ntheta
  x = np.zeros((N,4))
  x[:,0] = range(N)
  # Lay down left edges
  x[:,1] = sky_limits[0] +  x_box_size*np.mod(np.arange(N), nx) 
  # Lay down bottom edges
  x[:,2] = sky_limits[2] +  y_box_size*np.mod(np.floor(np.arange(N)/nx), ny)
  # Lay down orientation edges
  x[:,3] = 0.0 + theta_box_size * np.mod(np.floor(np.arange(N)/(nx*ny)), ntheta)
  # Add random offsets
  x[:,1] += np.random.uniform(0.,x_box_size, size=N)
  x[:,2] += np.random.uniform(0.,y_box_size, size=N)
  x[:,3] += np.random.uniform(0.,theta_box_size, size=N)
  return x

if __name__ == "__main__":
  number_passes = 3
  plt.figure(figsize=(15,15))
  plt.subplot(221)
  xA = generate_uniform_survey(number_passes)
  plot_survey(xA,'A')
  np.savetxt('A.txt',xA)
  
  xB = generate_uniform_survey(number_passes, rotate = True)
  plt.subplot(222)
  plot_survey(xB,'B')
  np.savetxt('B.txt',xB)
  
  xC = generate_uniform_survey(number_passes, offset = True)
  plt.subplot(223)
  plot_survey(xC,'C')
  np.savetxt('C.txt',xC)  
  
  xD = generate_random_survey(len(xA))
  plt.subplot(224)
  plot_survey(xD,'D')
  np.savetxt('D.txt',xD)
  plt.savefig("./surveys.png",bbox_inches='tight',pad_inches=0.)
  