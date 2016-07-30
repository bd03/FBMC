function [rx_output] = func_Analysis_Filter_Bank(y_sum,M,K,h,num_symbols,num_oqam_subsymbols,lp,noisy,amplitude)
%% func_Analysis_Filter_Bank
%
% Burak Dayi
%
% This function will implement FBMC RX Analysis Filter Bank 
%
% Created: 24-11-2014
%% Receiver
% the matrix that reshaped input samples would be stored in
receiver_input_1 = zeros(M,K+ceil(num_oqam_subsymbols/2)-1); 
receiver_input_2 = zeros(M,K+floor(num_oqam_subsymbols/2)-1);

% reshaping will be separated
% reshaping (joint implementation of delay chain & downsamplers)
% PPN1
for r=0:K+ceil(num_oqam_subsymbols/2)-1-1
    receiver_input_1(:,r+1) = (1/amplitude)*y_sum(1,1+r*M:M+r*M);
end

for r=0:K+floor(num_oqam_subsymbols/2)-1-1
    receiver_input_2(:,r+1) = (1/amplitude)*y_sum(1,1+M/2+r*M:M+M/2+r*M);
end


rx_poly_output_1 = zeros(M,K+ceil(num_oqam_subsymbols/2)-1+K-1); %output of PPN1 will be stored in this
rx_poly_output_2 = zeros(M,K+floor(num_oqam_subsymbols/2)-1+K-1); %output of PPN2 will be stored in this

jp_rx_poly_output = zeros(M,2*(K+num_symbols-1+K-1)); %output of polyphase filters will be stored in this
% jp_rx_fft_input = zeros(M,2*num_symbols);

% polyphase filter coefficients for archive
% ppb = zeros(M,K);

% polyphase filters are applied
for k=1:M
    b = h(M-k+1:M:lp+1); % related polyphase filter coefficients sieved
    % we use lp+1 b/c of the delay implemented in prototype filter design   
    rx_poly_output_1(k,:) = conv(receiver_input_1(k,:),b);
    rx_poly_output_2(k,:) = conv(receiver_input_2(k,:),b);
%     ppb(k,:) = b;
end

% rearrangement
len=K+ceil(num_oqam_subsymbols/2)-1+K-1+K+floor(num_oqam_subsymbols/2)-1+K-1;
rx_fft_input = zeros(M,len);
rx_fft_input(:,1:2:end)=rx_poly_output_1;
rx_fft_input(:,2:2:end)=rx_poly_output_2;

% save('ppb.mat','ppb');

% fft performed
rx_fft_output=fft(rx_fft_input)/sqrt(M);
% jp_rx_fft_output=fft(jp_rx_poly_output);

rx_output = zeros(M,len);
% jp_rx_output = zeros(M,2*(num_symbols+K-1+K-1));

% we convolve one sample with the entire filter. then at the receiver we
% convolve it with another filter. After contributions from other
% subcarrier branches were sieved out at the FFT due to orthogonality,
% we're only left with the sample simply convolved two times or in other
% words multiplied with (K+K-1) coefficients. The sumfactor is the sum of
% these coefficients that will enable us to normalize the sample that we
% get at the end. This factor is the same for all polyphase filter pairs.
% Therefore we are allowed to use any Aq-Bq pair.
% sumfactor = sum(conv(ppa(1,:),ppb(1,:)));
sumfactor = 0.6863;

for k=1:M
    beta = ((-1).^((k-1)*(0:(len-1))))*((-1).^((k-1)*K));

    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    
    rx_output(k,:) = rx_fft_output(k,:).*beta;
    % rx_output(k,:) = (1/sumfactor)*rx_fft_output(k,:).*beta;
end

% rx_output will be sent to subchannel processing block
% disp('+Receiver Block is processed.');