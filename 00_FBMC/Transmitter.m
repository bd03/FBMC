%% Transmitter
%
% Burak Dayi
%
% This script will perform transmitter block (synthesis filter bank).
%
% Dependencies: oqam_m - OQAM modulated message, num_symbols, K,M, lp
% Output: y - composite signal output
%
% Created: 02-03-2014

%disp('Transmitter Block')
y = zeros(1,K*M+(2*num_symbols-1)*M/2); %composite signal is in series. (k>=1)

% beta multiplier will be an exponential value multiplied by (-1)^kn where
% k is the subchannel number and n is subsymbol index. The subsymbol index
% we are using is 0 for the first element in a row and 1 is for the second
% element. Therefore it simplifies the general equation for the first
% element, we only need to compute exponential, in general. More
% simplifications also may be made.
beta = ones(1,2*num_symbols); % beta multiplier
ifft_input = ones(M,2*num_symbols); % ifft input


% polyphase filter & beta coefficients for archive
ppa = zeros(M,K);
betas = zeros(M,2*num_symbols);

for k=1:M
    bb = [exp(-j*2*pi*(k-1)*(lp-1)/(2*M)) ...
        ((-1)^(k-1))*exp(-j*2*pi*(k-1)*(lp-1)/(2*M))];
    
    if lp==K*M-1 %special treatment due to extra delay inserted 
        bb=[(-1)^((k-1)*K) (-1)^((k-1)*(K+1))];
    end
    
    beta =[];
    
    for c = 1:num_symbols
        beta = [beta bb];
    end
    
    betas(k,:)=beta;
    
    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    ifft_input(k,:)= oqam_m(k,:).*beta;
end

%ifft
tx_ifft_output=ifft(ifft_input);
tx_poly_output = zeros(M,2*(K+num_symbols-1)); %output of polyphase filters will be stored in this
%upsampler_output = zeros(M,M*(K+1)/2);



%polyphase filters are applied
for k=1:M
    a = h(k:M:lp+1); % related polyphase filter coefficients sieved
    % we use lp+1 b/c of the delay implemented in prototype filter design   
    %%%%%%%%%%tx_poly_output(k,:) = conv(a,tx_ifft_output(k,:));
    %we seperated the real and imag part chains as advised in primer
    %document
    tx_poly_output(k,:) = [conv(a,tx_ifft_output(k,1:2:2*num_symbols)) ...
        conv(a,tx_ifft_output(k,2:2:2*num_symbols))];
    ppa(k,:) = a;
    
end

save('ppa.mat','ppa');
save('betas.mat','betas');

% circular shift'i kodu guzellestirirken kullan simdi degil. simdi kolay ve
% acik yoldan yap.

%upsampling will prouduce
%upsampling, delay chain and summation are jointly implemented
% for r=0:1
%     y(1,1+r*M/2:K*M+r*M/2)= y(1,1+r*M/2:K*M+r*M/2)+tx_poly_output(:,r+1).';
% end


for r=0:K+num_symbols-1-1
    y(1,1+r*M:M+r*M) = y(1,1+r*M:M+r*M)+ tx_poly_output(:,r+1).';
    y(1,1+M/2+r*M:M+M/2+r*M) = y(1,1+M/2+r*M:M+M/2+r*M) + tx_poly_output(:,r+K+num_symbols).';
end

% y will be sent thru the channel