import ctypes


# import the C functions


# initialize the simulation parameters


# cal


# How to make an array from a python list
pyarr = [1, 2, 3, 4]
arr = (ctypes.c_int * len(pyarr))(*pyarr)