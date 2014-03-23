% test3.m
close all
clear all
clc

m = (randi(16,8,4)-1).'; %random samples generated
qam_m = qammod(m, 16, pi/2,'gray'); %built-in MATLAB qam modulation 
