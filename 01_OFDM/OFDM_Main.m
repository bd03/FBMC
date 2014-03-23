%% OFDM_Main
%
% Burak Dayi
%
% This is the main script that will call other scripts/functions

close all
clear all
clc
fprintf('--------------\n-----OFDM-----\n--------------\n\n');

%% Transmission
Config;
disp('+Configuration is obtained.');
Symbol_Creation;
disp('+Symbols are created.');
Transmitter;
disp('+Transmitter Block is processed.');

%% Channel
Channel;
disp('+Channel effects are applied.');

%% Receiver
Receiver;
disp('+Receiver Block is processed.');
Symbol_Estimation;
disp('+Symbol Estimation is done.');
Results;
disp('+Calculations and Statistics..');

