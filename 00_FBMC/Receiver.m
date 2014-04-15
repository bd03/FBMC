%% Receiver
%
% This will perform receiver block (analysis filter bank).
%
% Dependencies: y_ch - composite signal gone through channel response
% Output: rx_output - receiver block output signal
%
% Created: 02-03-2014

% disp('Receiver Block')

% The downsamplers following delay chain will form a Mx(K+1) matrix with
% channel output. Therefore, the delay chain and downsampler are jointly
% implemented as a vector reshape function.

% the matrix that reshaped input samples would be stored in
receiver_input = zeros(M,2*(K+num_symbols-1)); 

% reshaping (joint implementation of delay chain & downsamplers)
for r=0:K+num_symbols-1-1
%     y(1,1+r*M:M+r*M) = y(1,1+r*M:M+r*M)+ tx_poly_output(:,r+1).';
%     y(1,1+M/2+r*M:M+M/2+r*M) = y(1,1+M/2+r*M:M+M/2+r*M) + tx_poly_output(:,r+K+1).';
    receiver_input(:,r+1) = y_ch(1,1+r*M:M+r*M);
    receiver_input(:,r+K+num_symbols) = y_ch(1,1+M/2+r*M:M+M/2+r*M);    
end

rx_poly_output = zeros(M,2*(K+num_symbols-1+K-1)); %output of polyphase filters will be stored in this
rx_fft_input = zeros(M,2*num_symbols);

% polyphase filter coefficients for archive
ppb = zeros(M,K);

% polyphase filters are applied
for k=1:M
    b = h(M-k+1:M:lp+1); % related polyphase filter coefficients sieved
    % b = a(M-i+1);
    % we use lp+1 b/c of the delay implemented in prototype filter design   
    rx_poly_output(k,:) = [conv(receiver_input(k,1:(K+num_symbols-1)),b) ...
        conv(receiver_input(k,(K+num_symbols):2*(K+num_symbols-1)),b)];
%     rx_fft_input(k,:) = [rx_poly_output(k,K:K+num_symbols-1) ...
%         rx_poly_output(k, 2*K+num_symbols-2+K:2*K+num_symbols-2+K+num_symbols-1)];
    ppb(k,:) = b;
end

save('ppb.mat','ppb');

% fft performed
rx_fft_output=fft(rx_poly_output);

rx_output = zeros(M,2*(num_symbols+K-1+K-1));

% we convolve one sample with the entire filter. then at the receiver we
% convolve it with another filter. After contributions from other
% subcarrier branches were sieved out at the FFT due to orthogonality,
% we're only left with the sample simply convolved two times or in other
% words multiplied with (K+K-1) coefficients. The sumfactor is the sum of
% these coefficients that will enable us to normalize the sample that we
% get at the end. This factor is the same for all polyphase filter pairs.
% Therefore we are allowed to use any Aq-Bq pair.
sumfactor = sum(conv(ppa(1,:),ppb(1,:)));

for k=1:M
    bb = [exp(-j*2*pi*(k-1)*(lp+1)/(2*M)) ((-1)^(k-1))*exp(-j*2*pi*(k-1)*(lp+1)/(2*M))];
    
    if lp==K*M-1 %special treatment due to extra delay inserted 
        bb=[(-1)^((k-1)*K) (-1)^((k-1)*(K+1))];
    end
    
    beta =[];
    
    for c = 1:(num_symbols+K-1+K-1)
        beta = [beta bb];
    end
   
    
    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    
    
    rx_output(k,1:2:end-1) =  rx_fft_output(k,1:(num_symbols+K-1+K-1));
    rx_output(k,2:2:end) =  rx_fft_output(k,(num_symbols+K-1+K-1)+1:2*(num_symbols+K-1+K-1));
    
    rx_output(k,:) = rx_output(k,:).*beta;
    
%     rx_output(k,:) = [];(1/sumfactor)*[sum(rx_fft_output(k,1:2*K-1)) ...
%         sum(rx_fft_output(k,2*K:2*(2*K-1)))].*beta; %(1/sumfactor)*
    
end

% rx_output will be sent to subchannel processing block