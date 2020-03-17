% Script for plotting the results from auto-coded ilqr

clear;
close all;
clc;


% Read the data from the csv
c_data = readmatrix("c_ilqr_run1.csv");
x = c_data(:,2);
y = c_data(:,3);
theta = c_data(:,4);
v = c_data(:,5);

figure(1)
plot(x, y);