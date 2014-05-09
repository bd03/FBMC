%% Symbol_Creation
%
% Burak Dayi
%
% This script will create the sequence of sufficient bit and then constructs
% FBMC symbols with the modulation scheme defined in configuration file.
%
% Dependencies: modulation, num_symbols, bits_per_sample, M, num_bits
% Output: qam_m - QAM modulated message, bits - transmitted bits
%
% Created: 02-03-2014

% disp('Symbol Creation')

% bit sequence creation
bits = randi(2,1,num_bits)-1;
% load('bits.mat');
%save('bits.mat','bits');

m = reshape(bi2de(reshape(bits,bits_per_sample,M*num_symbols).','left-msb'),M,num_symbols);
%m = (randi(modulation,num_symbols,num_samples)-1).'; %random samples generated

qam_m = qammod(m, modulation, pi/2,'gray'); %built-in MATLAB qam modulation 
