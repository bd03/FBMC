%% Symbol_Creation
%
% Burak Dayi
%
% This script will create the sequence of sufficient bit and then constructs
% OFDM symbols with the modulation scheme defined in configuration file.
%

% normalization constants
normalization = [2 10 42 170];
q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.

% bit sequence creation
bits = randi(2,1,num_bits)-1;
upperside = zeros(indices(1)-1,num_symbols);
lowerside = zeros(M-indices(2),num_symbols);
m_data = reshape(bi2de(reshape(bits,bits_per_sample,allocated_subchannels*num_symbols).','left-msb'),allocated_subchannels,num_symbols);
m = [upperside; m_data; lowerside];
%m = (randi(modulation,num_symbols,allocated_subchannels)-1).'; %random samples generated
qam_m = qammod(m, modulation, pi/2,'gray'); %built-in MATLAB qam modulation 
qam_m = qam_m/sqrt(normalization(find(modulation==q_arr)));
% qam_m =pskmod(m,modulation,pi/2,'gray');