function [bits_est,m_est] = func_Symbol_Estimation(num_bits,bits_per_sample,N,rx_output,modulation,amp)
	%% Symbol_Estimation
%
% Burak Dayi
%
% Created on: 05-02-2015

% normalization constants
normalization = [2 10 42 170];
num_symbols = size(rx_output,2);

m_est = zeros(N,num_symbols);

% Estimation of the symbols
for k=1:N
    m_est(k,:) = qamdemod(sqrt(normalization(log2(modulation)/2))*rx_output(k,:)/amp,modulation,pi/2,'gray');
%     m_est(k,:) = pskdemod(rx_output(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
bits_est=reshape(de2bi(reshape(m_est,num_symbols*N,1),'left-msb').',1,num_bits);