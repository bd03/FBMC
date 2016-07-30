%% SCFDMA_Main
%
% Burak Dayi
%
% This is the main script that will call other scripts/functions
%
% Created on 27-01-2015

close all
clear all
clc
fprintf('---------------\n----SC-FDMA----\n---------------\n\n');

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