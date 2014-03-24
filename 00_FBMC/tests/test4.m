%test4
%% pskmod test

close all
clear all
clc

m = (randi(16,8,4)-1).'; %random samples generated
psk_m = pskmod(m, 16, pi/2,'gray'); %built-in MATLAB qam modulation

psk_demod = pskdemod(psk_m,16,pi/2,'gray');

err=m-psk_demod;
