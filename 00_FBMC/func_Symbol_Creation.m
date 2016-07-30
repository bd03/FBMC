function [m,bits,output] = func_Symbol_Creation(num_bits,bits_per_sample, len, num_symbols,modulation)
%% func_Symbol_Creation
%
% Burak Dayi
%
% This function will create the sequence of sufficient bit and then constructs
% FBMC symbols with the modulation scheme defined in configuration file.
%
% Dependencies: modulation, num_symbols, bits_per_sample, len, num_bits
% Output: output - QAM modulated message
%
% Created: 21-11-2014

% normalization constants
normalization = [2 10 42 170];
q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.

% bit sequence creation
bits = randi(2,1,num_bits)-1;
% load('bits.mat');
%save('bits.mat','bits');

m = reshape(bi2de(reshape(bits,bits_per_sample,len*num_symbols).','left-msb'),len,num_symbols);
%m = (randi(modulation,num_symbols,num_samples)-1).'; %random samples generated

qam_m = qammod(m, modulation, pi/2,'gray'); %built-in MATLAB qam modulation 
output = qam_m/sqrt(normalization(find(modulation==q_arr)));