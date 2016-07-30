function [A,yy] = func_Synthesis_Filter_Bank(oqam_m,M,K,h,num_oqam_subsymbols,lp,noisy,amp)
%% func_Synthesis_Filter_Bank.m
%
% Burak Dayi
%
% This function will implement FBMC Synthesis Filter Bank 
%
% Created: 24-11-2014

%% Transmitter
%disp('Transmitter Block')

y = zeros(1,K*M+(num_oqam_subsymbols-1)*M/2); %composite signal is in series. (k>=1)
jp_y = zeros(1,K*M+(num_oqam_subsymbols-1)*M/2); 

% beta multiplier will be an exponential value multiplied by (-1)^kn where
% k is the subchannel number and n is subsymbol index. The subsymbol index
% we are using is 0 for the first element in a row and 1 is for the second
% element. Therefore it simplifies the general equation for the first
% element, we only need to compute exponential, in general. More
% simplifications also may be made.
beta = ones(1,num_oqam_subsymbols); % beta multiplier
ifft_input = ones(M,num_oqam_subsymbols); % ifft input


% polyphase filter & beta coefficients for archive
% ppa = zeros(M,K);
betas = zeros(M,num_oqam_subsymbols);

for k=1:M
    
    beta = ((-1).^((k-1)*(0:(num_oqam_subsymbols-1))))*((-1).^((k-1)*K)); %special treatment for prototype filters of length KM-1
    
    betas(k,:)=beta;
    
    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    ifft_input(k,:)= oqam_m(k,:).*beta;
end

%ifft
tx_ifft_output=ifft(ifft_input)*sqrt(M);

%polyphase filters are applied
for k=1:M
    a = h(k:M:lp+1); % related polyphase filter coefficients sieved
    % we use lp+1 b/c of the delay implemented in prototype filter design   
    %we seperated the real and imag part chains as advised in primer
    %document
    tx_poly_output_1(k,:) = conv(a,tx_ifft_output(k,1:2:num_oqam_subsymbols)); %output of the 1st PPN
    tx_poly_output_2(k,:) =conv(a,tx_ifft_output(k,2:2:num_oqam_subsymbols)); %output of the 2nd PPN
    ppa(k,:) = a;
end

save('ppa.mat','ppa');
save('betas.mat','betas');

% contribution from PPN1
for r=0:K+ceil(num_oqam_subsymbols/2)-1-1
    y(1,1+r*M:M+r*M) = y(1,1+r*M:M+r*M)+ tx_poly_output_1(:,r+1).';
end

% contribution from PPN2
for r=0:K+floor(num_oqam_subsymbols/2)-1-1
    y(1,1+M/2+r*M:M+M/2+r*M) = y(1,1+M/2+r*M:M+M/2+r*M) + tx_poly_output_2(:,r+1).';
end
yy=y;
A= NaN;
if noisy
    % current output power:
    e_pow = sum(abs(y.*y));
    % e_pow=sum(abs(y));
    % amp = sqrt(y_pow_req/y_pow) %% however, one can also use absolute power
    % amp = amp/y_pow
    A = amp/sqrt(e_pow);
    yy = A*y;
end



%y = user_amps*y;

% y will be sent thru the channel