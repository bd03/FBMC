%% Symbol_Creation
%
% Burak Dayi
%
% This script will create the sequence of sufficient bit and then constructs
% OFDM symbols with the modulation scheme defined in configuration file.
%
% Created on 26-01-2015


% bit sequence creation
bits = randi(2,1,num_bits)-1;
m = reshape(bi2de(reshape(bits,bits_per_sample,N*num_symbols).','left-msb'),N,num_symbols); % for now the array size is N*num_symbols
%m = (randi(modulation,num_symbols,num_samples)-1).'; %random samples generated
qam_m = qammod(m, modulation, pi/2,'gray')/sqrt(normalization(log2(modulation)/2)); %built-in MATLAB qam qam 
% qam_m =pskmod(m,modulation,pi/2,'gray');