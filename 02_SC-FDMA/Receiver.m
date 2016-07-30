%% Receiver
% Burak Dayi
%
% This script will perform receiver block.
%
% Created on: 30-01-2015

%% Cyclic Prefix Removal
rx_serial = zeros(1,M*num_symbols);

for i=1:num_symbols
    rx_serial(1,1+(i-1)*M:i*M)= y_ch(1+cp_length+(i-1)*(M+cp_length):i*(M+cp_length));
end

%% Serial to Parallel Conversion

rx_fft_input = reshape(rx_serial,M,num_symbols);


%% FFT
rx_fft_output = fft(rx_fft_input)/sqrt(M);%/sqrt(N/(M^2));%/sqrt(N);
rx_demapped = subcarrier_demapping(rx_fft_output,Q,N,start,mapping)/sqrt(M/N); % Subcarrier demapping
rx_output = ifft(rx_demapped)*sqrt(N);%/(N/M); % N-point DFT first.