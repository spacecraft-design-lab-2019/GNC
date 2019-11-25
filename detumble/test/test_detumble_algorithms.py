'''
Script for testing detumble algorithms
''' 
import detumble_cpp as dcpp
import numpy as np
import sys
import pytest
from random import gauss
#from GNC.cmake_build_debug import detumble_algorithms_cpp as dtcpp


def make_rand_vector(dims):
    '''
    :param dims: dimension of randomly distributed vector
    :return: a vector in 3 space
    '''
    vec = [gauss(0, 1) for i in range(dims)]
    mag = sum(x**2 for x in vec) ** .5
    return [x/mag for x in vec]


def test_bias_estimation():
    # known bias
    bias_true = np.array([40000, 20000, 10000])

    # generate fake biased measurements
    m = 200
    B_mat = np.zeros((m, 3))
    mu, sigma = 45000., 5000.  # mean and standard deviation
    B_mags = np.random.normal(mu, sigma, m)

    for i in range(m):
        rand_vec = np.array(make_rand_vector(3))
        B_true = np.transpose(bias_true) + B_mags[i]*rand_vec
        B_mat[i,:] = np.transpose(B_true)
    bias_est = dcpp.get_bias_estimate(B_mat)

    abs_error = (bias_est-bias_true)
    rel_error = np.linalg.norm(abs_error)/np.linalg.norm(bias_true)


    np.testing.assert_array_less(rel_error,.05)


    
    
    