%% Symbol_Estimation
%
% This will recombine the QAM samples and extract estimated symbol from
% those samples and transmitted bit estimation from those symbols.
%
% Dependencies: oqam_demod - OQAM demodulated signal samples, M,
% num_symbols, bits_per_sample, num_bits
% Output: m_est - estimated symbols, bits_est - estimated bits
%
% Created: 17-03-2014

% disp('Symbol Estimation')

qam_est = zeros(M,num_symbols);
m_est = zeros(M,num_symbols);

% Estimation of the symbols
for k=1:M
    qam_est(k,:)= (1/sumfactor)*(oqam_demod(k,1:2:2*num_symbols-1)+oqam_demod(k,2:2:2*num_symbols));
    m_est(k,:) = qamdemod(qam_est(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
%bits_est=de2bi(reshape(m_est,1,num_symbols*M));
bits_est=reshape(de2bi(reshape(m_est,num_symbols*M,1),'left-msb').',1,num_bits);