%% Transmitter
%
% Burak Dayi
%
% This script will perform transmitter block.
%
%
% Last updated: 18-03-2014

tx_ifft_output = ifft(qam_m);

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

%tx_output is transmitted through channel