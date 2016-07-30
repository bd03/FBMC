%% Symbol_Estimation
%
% Burak Dayi
%
% Created on: 30-01-2015

m_est = zeros(N,num_symbols);

% Estimation of the symbols
for k=1:N
    m_est(k,:) = qamdemod(sqrt(normalization(log2(modulation)/2))*rx_output(k,:),modulation,pi/2,'gray');
%     m_est(k,:) = pskdemod(rx_output(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
bits_est=reshape(de2bi(reshape(m_est,num_symbols*N,1),'left-msb').',1,num_bits);