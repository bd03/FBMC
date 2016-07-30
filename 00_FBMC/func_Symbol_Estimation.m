function [bits_est,m_est] = func_Symbol_Estimation(oqam_demod,num_bits,allocated_subchannels,num_symbols,modulation)
%% func_Symbol_Estimation
%
% Burak Dayi
%
% This function will implement conversion of QAM symbols into bitstream
%
% Created: 24-11-2014

%% Symbol_Estimation
%
sumfactor = 0.6863;
% normalization constants
normalization = [2 10 42 170];
q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.
% num_samples = allocated_subchannels*num_symbols; %number of samples in a vector
% bits_per_sample = log2(modulation); %num of bits carried by one sample
% num_bits = num_samples*bits_per_sample; % total number of bits transmitted

qam_est = zeros(allocated_subchannels,num_symbols);
m_est = zeros(allocated_subchannels,num_symbols);

% Estimation of the symbols
for k=1:allocated_subchannels
    qam_est(k,:)= (1/(sumfactor))*(oqam_demod(k,1:2:2*num_symbols-1)+oqam_demod(k,2:2:2*num_symbols))*sqrt(normalization(find(modulation==q_arr)));
    % qam_est(k,:)= (1/(sumfactor*amplitude))*(oqam_demod(k,1:2:2*num_symbols-1)+oqam_demod(k,2:2:2*num_symbols))*sqrt(normalization(find(modulation==q_arr)));
    m_est(k,:) = qamdemod(qam_est(k,:),modulation,pi/2,'gray');
end

% Estimation of transmitted bits
%bits_est=de2bi(reshape(m_est,1,num_symbols*M));
bits_est=reshape(de2bi(reshape(m_est,num_symbols*allocated_subchannels,1),'left-msb').',1,num_bits);
% disp('+Symbol Estimation is done.');