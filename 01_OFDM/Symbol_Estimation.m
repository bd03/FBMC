%% Symbol_Estimation
%
% Burak Dayi
%
%

m_est = zeros(M,num_symbols);

% Estimation of the symbols
for k=1:M
    m_est(k,:) = qamdemod(rx_output(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
%bits_est=de2bi(reshape(m_est,1,num_symbols*M));
bits_est=reshape(de2bi(reshape(m_est,num_symbols*M,1),'left-msb').',1,num_bits);