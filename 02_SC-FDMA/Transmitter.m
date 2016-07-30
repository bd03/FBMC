%% Transmitter
%
% Burak Dayi
%
% This script will perform transmitter block.
%
%
% Created on 26-01-2015

tx_dft_output = fft(qam_m)/sqrt(N); % N-point DFT first.
tx_mapped = subcarrier_mapping(tx_dft_output,Q,M,start,mapping); % Subcarrier mapping
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

% tx_output is normalized
tx_output = sqrt(noise_floor_wattsHz * 10^(0.1*EbN0)*num_bits)*tx_output;

%tx_output is transmitted through channel