function [y,num_oqam_subsymbols] = func_FBMC_Transmitter(M,K,h,num_symbols,low_boundary,high_boundary,qam_m,lp,preamble_enabled,preamble,syms_per_frame)
%% func_FBMC_Transmitter
%
% Burak Dayi
%
% This function will implement FBMC TX starting from OQAM preprocessing 
%
% Created: 21-11-2014

%% OQAM Preprocessing
% disp('OQAM preprocessing')

jn = (j.^((1:2*num_symbols)-1)).'; % j to the power of n
%for odd subchannels flipped version of this will be used.

oqam_m = zeros(M, 2*num_symbols); %the matrix that will store modulated signal

symbol_matrix_index = 1; % keeps the index number in qam_m matrix
for k=low_boundary:high_boundary
    theta = jn*(j^(k-1)); %theta multiplier
    if ~mod((k-1),2) % even k --> index spans from 1 to 256
        real_parts = upsample(real(qam_m(symbol_matrix_index,:).'),2);
        imag_parts = circshift(upsample(imag(qam_m(symbol_matrix_index,:).'),2),1);        
    else % odd k
        real_parts = circshift(upsample(real(qam_m(symbol_matrix_index,:)).',2),1);
        imag_parts = upsample(imag(qam_m(symbol_matrix_index,:).'),2);           
    end
    oqam_m(k,:) = (real_parts+imag_parts).*theta;
    symbol_matrix_index = symbol_matrix_index +1;
end

jp = zeros(size(oqam_m));

if preamble_enabled
    % preamble insertion
    for i=num_frames:-1:1
        oqam_m = [oqam_m(:,1:2*(i-1)*syms_per_frame) preamble oqam_m(:,2*(i-1)*syms_per_frame+1:end)];
        jp = [jp(:,1:2*(i-1)*syms_per_frame) preamble jp(:,2*(i-1)*syms_per_frame+1:end)];
    end
end
% 
% %%!!!!!!!!!!!!!!!!hack!!!!!!!!!!!!!!!!!!!!!!!!!!
%num_symbols = num_symbols+num_frames*2; %num_frames= no of preambles
% num_oqam_subsymbols = 2*num_symbols+3*num_frames; %num_frames= no of preambles
num_oqam_subsymbols = size(oqam_m,2);
%in this file all the computation will depend on num_oqam_subsymbols.

% %soru: is preamble oqam modulated or to be modulated???????????????
% %ans: the preamble takes 3 subsymbol duration with given setup. that
% %indicates that the preamble will be appended to the oqam modulated data



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
    
    beta = ((-1).^((k-1)*(1:num_oqam_subsymbols)))*((-1).^((k-1)*K)); %special treatment for prototype filters of length KM-1
    
    betas(k,:)=beta;
    
    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    ifft_input(k,:)= oqam_m(k,:).*beta;
end

%ifft
tx_ifft_output=ifft(ifft_input);

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

% y will be sent thru the channel