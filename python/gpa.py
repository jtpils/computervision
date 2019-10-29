#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
    Last reviewed: 26/10/2019
    author: G. A. Moreira
'''

import scipy
import numpy as np
import matplotlib.pyplot as plt

from math import radians

from mpl_toolkits import mplot3d
from sklearn.neighbors import NearestNeighbors


def rotmat(alpha, beta, gamma):
    '''
        Builds a numpy.ndarray rotation matrix with the provided rotation
        angles, (in degrees). Gamma is the angle about the x-axis, beta about
        the y-axis and alpha about the z-axis.
    '''
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)

    Rx = np.array([[1, 0, 0],
                   [0, np.cos(gamma), -np.sin(gamma)],
                   [0, np.sin(gamma), np.cos(gamma)]])

    Ry = np.array([[np.cos(beta), 0, np.sin(beta)],
                   [0, 1, 0],
                   [-np.sin(beta), 0, np.cos(beta)]])

    Rz = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                   [np.sin(alpha), np.cos(alpha), 0],
                   [0, 0, 1]])

    R = np.matmul(Rz,Ry,Rx)
    
    return R



def mutual_nearest(clouds):
    '''
        Receives any number (N) of point clouds as an input. The point clouds 
        should be type numpy.ndarray with the first dimension being the points
        (P) and the second the coordinates (x,y,z) e.g.
        Outputs a numpy.ndarray of dimensions (Q,N). Q is the number of groups
        of mutually nearest neighbours and N are all the point clouds. This
        array contains the indices to select the points from each cloud. So 
        the 4th row of the array contains N indices to select the points from
        each cloud that form a mutually nearest neighbor group.
    '''    
    n_clouds = len(clouds)
    points_per_cloud = [cloud.shape[0] for cloud in clouds]
    n_points = np.sum(points_per_cloud)    
    correspondences = np.zeros([n_points, n_points])
    corr = np.zeros([1000, n_clouds, n_clouds], dtype='uint16')
    
    for i in range(n_clouds):
        for j in range(n_clouds):
            nbrs = NearestNeighbors(1, algorithm='ball_tree').fit(clouds[j])                                       
            _, indices_j = nbrs.kneighbors(clouds[i])
            corr[:, i, j] = indices_j.flatten()
        
    for i in range(n_clouds):
        for j in range(n_clouds):
            match = np.where(np.arange(0,corr.shape[0]) - corr[corr[:,i,j],j,i] == 0)
            index_i = corr[match,i,i].flatten() + i*points_per_cloud[i]
            index_j = corr[match,i,j].flatten() + j*points_per_cloud[j]
            correspondences[(index_i, index_j)] = 1
        
    return correspondences


def procrustes(A, B):
    '''
        Solves the Procrustes problem using the closed-form solution presented
        in Schonemann (1977). Given two (p,q) matrices A, B it finds the scalar
        S, the rotation matrix R and the translation vector T that minimize the
        error E = s*A*R + j*T' - B (i.e. the Frobenius norm of E). The solution
        involves the use of singular value decomposition (SVD) and is fully
        detailed in the article.
    '''
    
    centroid_B = np.mean(B, axis=0)
    centroid_A = np.mean(A, axis=0)
    
    A_centered = A - centroid_A
    B_centered = B - centroid_B
    
    C = np.matmul(A_centered.T, B_centered)
    
    V, _, W = np.linalg.svd(C)
    
    R = np.matmul(V, W).T
    detR = np.linalg.det(R)
    
    if detR < 0:
        reflex = np.diag([1,1,detR])
        R = np.matmul(np.matmul(V, reflex), W)
        
    T = centroid_B - np.matmul(R, centroid_A)
    s = 1
    
    return s, R, T  



def cloudviz(clouds, **kwargs):
    '''
        Receives a list of numpy.ndarray point clouds and represents them in a
        3D coordinate system.
    '''
    plt.figure(facecolor='0')
    ax = plt.axes(projection='3d')
    ax.set_facecolor('black')
    for i in range(len(clouds)):
        ax.plot(clouds[i][:,0],clouds[i][:,1],clouds[i][:,2], 'o', markersize=1)
    ax.legend(['Point cloud ' + str(i) for i in range(len(clouds))])
    plt.autoscale('off')
    plt.axis('off')
    plt.grid(b=None)
    plt.show()
    


def SRTtransform(cloud, s, R, T):
    '''
        Transforms a point cloud provided as a numpy.ndarray, according to a
        scale factor s, a rotation R and a translation T. The cloud should
        have dimensions (p,q) where p is the number of points and q is the 
        number of dimensions. s is a float, R a (3,3) numpy.ndarray and T is
        a (1,3) numpy.ndarray.
    '''
    rotated = np.matmul(R, cloud.T).T 
    scaled = s*rotated
    translated = scaled + T.reshape([1,3])

    return translated


class Centroid:
    '''
    '''
    def __init__(self, groups, stacked_clouds, points_per_cloud):
        '''
        '''
        self.centroid = []
        self.points_per_cloud = points_per_cloud
        self.groups = groups
        for i in range(len(self.groups)):
             points = stacked_clouds[self.groups[i]]
             mean = np.mean(points, axis=0)
             self.centroid.append(mean)
        
        self.centroid = np.array(self.centroid)


    def filterByCloud(self, cloud_number):
        '''
        '''
        start_index = np.sum(self.points_per_cloud[:cloud_number])
        end_index = np.sum(self.points_per_cloud[:1+cloud_number])
        cloud_points = []
        centroid_points = []
        for i in range(len(self.groups)):
            group = self.groups[i]
            condition = np.logical_and(group >= start_index, group < end_index)
            condition_true = np.where(condition)[0]
            if len(condition_true) == 1:
                centroid_points.append(self.centroid[i,:])
                cloud_points.append(stacked_clouds[group[condition_true[0]]])

        return np.array(cloud_points), np.array(centroid_points)
    
    
#%%


mat = scipy.io.loadmat('bunny.mat')
clouds_ = [mat['modelSimply'][0][i][0] for i in range(3)]



#%%
# Register the point clouds using Procrustes within a framework similar to ICP

clouds = clouds_.copy()

MAX_ITER = 50
N_CLOUDS = len(clouds)
POINTS_PER_CLOUD = [cloud.shape[0] for cloud in clouds]
    
for i in range(MAX_ITER):
    print('Iteration: ' + str(i))
    
    stacked_clouds = np.vstack(clouds)
    correspondence_map = mutual_nearest(clouds)
    vals, idx_start, count = np.unique(correspondence_map, axis=0, return_counts=True, return_index=True)
    points_per_group = np.sum(vals, axis=1)
    correspondence_map_unique = vals[np.logical_and(points_per_group==count, count>1)]
    
    group_id, point_indices = np.array(np.where(correspondence_map_unique==1), dtype='uint16')
    split_indices = np.unique(group_id, return_index=True)[1]
    groups = np.split(point_indices, split_indices)

    err_before = 0
    err_after = 0
    
    for k in range(N_CLOUDS):
        print('   Point cloud: ' + str(k))
        centroid = Centroid(groups, stacked_clouds, POINTS_PER_CLOUD)
        A, B = centroid.filterByCloud(k)
        print('      No. Points used (centroid): ' + str(A.shape[0]))
        
        cloud_error_before = np.sum(np.square(A-B))
        err_before += cloud_error_before
        print('      Cloud-Centroid MSE (before): ' + str(round(err_before/A.shape[0],3)))
        
        s, R, T = procrustes(A,B)
        clouds[k] = SRTtransform(clouds[k],s,R,T)

        cloud_error_after = np.sum(np.sum(np.square(SRTtransform(A,s,R,T)-B))) 
        err_after += cloud_error_after 
        print('      Cloud-Centroid MSE (after): ' + str(round(err_after/A.shape[0],3)))
        
    err_before /= stacked_clouds.shape[0]
    err_after /= stacked_clouds.shape[0]
    
    print('MSE before-after: ' + str(round(err_before,4)) +\
           str('  ' ) + str(round(err_after,4)))




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
