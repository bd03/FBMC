%% Symbol_Estimation
%
% Burak Dayi
%
%

% normalization constants
normalization = [2 10 42 170];
q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.

m_est = zeros(allocated_subchannels,num_symbols);

% norm_rx_output=rx_output(k,:)*sqrt(normalization(find(modulation==q_arr)));

% Estimation of the symbols
for k=indices(1):indices(2)
	norm_rx_output(k,:)= rx_output(k,:)*sqrt(normalization(find(modulation==q_arr)))/amplitude;
    m_est(k-indices(1)+1,:) = qamdemod(norm_rx_output(k,:),modulation,pi/2,'gray');
%     m_est(k,:) = pskdemod(rx_output(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
%bits_est=de2bi(reshape(m_est,1,num_symbols*M));
bits_est=reshape(de2bi(reshape(m_est,num_symbols*allocated_subchannels,1),'left-msb').',1,num_bits);