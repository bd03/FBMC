%% FBMC_Main
%
% Burak Dayi
%
% This is the main script that will call other scripts/functions.
%
% This will implement Filterbank Multicarrier transmission and reception
% with modulation and demodulations steps described in PHYDYAS documents
% Deliverable 5.1 and FBMC Physical Layer: A primer.
%
% Notes:
% Configuration will allow only subcarrier sizes of a power of two so as to
% be able to use FFT in the implementation
%
% The prototype filter used here is the PHYDYAS filter described in
% Deliverable 5.1: Prototype filter and structure optimization
% The implementation of the blocks follows what is described in Deliverable
% 5.1 transmultiplexer architecture.
%
% Channel implementation
% Subchannel processing
% Calculations and statistics
%
% Created: 02-03-2014

close all
clear all
clc
fprintf('--------------\n-----FBMC-----\n--------------\n\n');

%% Transmission 
Config;
disp('+Configuration is obtained.');
Prototype_filter;
disp('+Prototype filter is designed.');
Symbol_Creation;
disp('+Symbols are created.');
OQAM_Preprocessing;
disp('+OQAM Preprocessing is done.');
Transmitter;
disp('+Transmitter Block is processed.');

%% Channel
Channel;
disp('+Channel effects are applied.');

%% Reception
Receiver;
disp('+Receiver Block is processed.');
Subchannel_processing;
disp('+Subchannel processing is done.');
OQAM_Postprocessing;
disp('+OQAM Preprocessing is done.');
Symbol_Estimation;
disp('+Symbol Estimation is done.');
Results;
disp('+Calculations and Statistics..');