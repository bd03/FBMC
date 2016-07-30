function [yy,A] = func_Transmitter(N,M,Q,mapping,start_index,qam_m,cp_length,noisy,amp)
%% func_Transmitter.m
%
% Burak Dayi
%
% This function will implement SC-FDMA Transmitter
%
% Created: 05-02-2015

num_symbols = size(qam_m,2);

tx_dft_output = fft(qam_m)/sqrt(N); % N-point DFT first.
tx_mapped = subcarrier_mapping(tx_dft_output,Q,M,start_index,mapping); % Subcarrier mapping
tx_ifft_output = ifft(tx_mapped)/sqrt(N/(M^2)); % M-point IFFT


%% Parallel to Serial Conversion
tx_serial = reshape(tx_ifft_output,1,M*num_symbols);

%% Cyclic prefix addition
tx_output = zeros(1,M*(num_symbols+cp_length));

for i=1:num_symbols
    tx_output(1,1+(i-1)*(M+cp_length):(i-1)*(M+cp_length)+cp_length) = ...
        tx_serial(1,i*M-cp_length+1:i*M); % prefix extension
    tx_output(1,1+(i-1)*(M+cp_length)+cp_length:(i-1)*(M+cp_length)+cp_length+M) = ...
        tx_serial(1,1+(i-1)*M:i*M); % the symbol
end

yy = tx_output;
% out2 = tx_output;
A= NaN;
if noisy
    % current output power:
    e_pow = sum(abs(tx_output.*tx_output));
    % e_pow=sum(abs(y));
    % amp = sqrt(y_pow_req/y_pow) %% however, one can also use absolute power
    % amp = amp/y_pow
    A = amp/sqrt(e_pow);
    yy = A*tx_output;
end