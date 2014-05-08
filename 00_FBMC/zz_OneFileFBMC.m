%% zz_OneFileFBMC.m
%
% Burak Dayi

% Sending one frame

%
% Created: 17-03-2014

close all
clear all
clc
fprintf('--------------\n-----FBMC-----\n--------------\n\n');

%% Transmission 
%% Config

% 2.b Main mode parameters
%---- General filterbank parameters ----%
K = 4; % overlapping factor 
M = 512; % number of subcarriers
% num_frames = 0; % number of data frames in each FBMC block
num_symbols = 50; % number of symbols sent back to back in one transmission
num_samples = M; %number of samples in a vector
modulation = 4; %4-, 16-, 64-, 128-, 256-QAM
bits_per_sample = log2(modulation); %num of bits carried by one sample
num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
lp = K*M-1; % filter length
delay = K*M+1-lp; %delay requirement

%---- Prototype filter frequency coefficients----%
%---------------------------------
% This part should not be altered!
%---------------------------------
% K will select the row, the last column is background noise power in dB 
% (that might be come in handy in future) for reference
P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];

%---- Preamble creation IAM method----%
preamble = [repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];

%---- Channel settings ----%
% noise settings
noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
SNR = 10; % SNR of the channel. ideal=0 to see the effects on channel

%rayleigh channel settings
fading = 0; % set 0 for distortionless channel, set 1 for rayleigh channel
bw = 5e+6; % Transmission Bandwidth
channel_profiles = ['EPA' 'EVA' 'ETU']; % Valid channel profile selections
profile ='EPA'; %Channel profile
[delay_a pow_a] = LTE_channels (profile,bw);
ch_resp = rayleighchan(1/bw,10,delay_a,pow_a); %channel model
ch_h.storeHistory = 1;
ch_h.storePathGains =1;

%---- Equalizer settings ----%
eq_select = 4; % selection of equalizer type 1: one tap, 
% 2: three tap w/ geometric interp, 3: three tap w/ linear interp
% 4: no equalizer


%---- Parameter check ---%
if K>4 || K<2
K
    error('Range of integer K is [2 4]');
end

if ~~mod(log2(M),1)
    M
    error('The number of subcarriers M should be a power of 2.');
end

if ~(lp == K*M || lp == K*M-1 || lp == K*M+1)
    lp
    error('Only filter sizes of KM, KM+1, KM-1 are suported.');
end

if ~~mod(log2(modulation),1)
    modulation
    error('Only modulation=2^m-QAM schemes are supported.');
end

if eq_select>4 || eq_select<1 || mod(eq_select,1)~=0
    eq_select
    error('eq_select should be an integer in range [1 4].');
end


% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if noisy
    disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d,\nfading=%d, equalizer=%d, estimation=''POP'', SNR=%d dB', K,M,num_symbols,modulation,num_bits,fading,eq_select,SNR));
else
    disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d,\nfading=%d, equalizer=%d, estimation=''POP'', SNR=Ideal', K,M,num_symbols,modulation,num_bits,fading,eq_select));
end
disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
%pause;

disp('+Configuration is obtained.');
Prototype_filter;
disp('+Prototype filter is designed.');
%% Symbol_Creation

% bit sequence creation
bits = randi(2,1,num_bits)-1;

m = reshape(bi2de(reshape(bits,bits_per_sample,M*num_symbols).','left-msb'),M,num_symbols);
%m = (randi(modulation,num_symbols,num_samples)-1).'; %random samples generated
qam_m = qammod(m, modulation, pi/2,'gray'); %built-in MATLAB qam modulation 

disp('+Symbols are created.');
%% OQAM_Preprocessing

jn = (j.^((1:2*num_symbols)-1)).'; % j to the power of n
%for odd subchannels flipped version of this will be used.

oqam_m = zeros(M, 2*num_symbols); %the matrix that will store modulated signal

for k=1:M
    theta = jn*(j^(k-1)); %theta multiplier
    if ~mod((k-1),2) % even k
        real_parts = upsample(real(qam_m(k,:).'),2);
        imag_parts = circshift(upsample(imag(qam_m(k,:).'),2),1);        
    else % odd k
        real_parts = circshift(upsample(real(qam_m(k,:)).',2),1);
        imag_parts = upsample(imag(qam_m(k,:).'),2);           
    end
    oqam_m(k,:) = (real_parts+imag_parts).*theta; 
end

% preamble insertion
oqam_m = [preamble oqam_m];
% 
% %%!!!!!!!!!!!!!!!!hack!!!!!!!!!!!!!!!!!!!!!!!!!!
num_symbols = num_symbols+1;
% %soru: is preamble oqam modulated or to be modulated???????????????
% %ans: the preamble takes 3 subsymbol duration with given setup. that
% %indicates that the preamble will be appended to the oqam modulated data
% %diger cvp: D3.1'de oyle demiyor ama


disp('+OQAM Preprocessing is done.');
%% Transmitter

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
disp('+Transmitter Block is processed.');

%% Channel

if fading
    y_filtered = filter(ch_resp,y);
else
    y_filtered = y;%filter(1,1,y);
end

if ~noisy
    y_ch = y_filtered;
else
    y_ch = awgn(y_filtered,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
end
disp('+Channel effects are applied.');

%% Reception
%% Receiver

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
disp('+Receiver Block is processed.');

%% Subchannel_processing

% remove preamble
scp_preamble = rx_output(:,2*K-1:2*K);
scp_data = rx_output;

% %%!!!!!!!!!!!!!!!!hack!!!!!!!!!!!!!!!!!!!!!!!!!!
num_symbols = num_symbols-1;
% 
% %estimation of channel
ch_resp_est = scp_preamble(:,1)./(preamble(:,1)*sumfactor);

if eq_select == 1
    % one tap equalizer
    for i=2*K+1:2*K+2*num_symbols
        sp_output(:,i-2*K) = scp_data(:,i)./ch_resp_est;
    end
elseif eq_select == 2 || eq_select == 3
    % % three tap equalizer
    % % coef computation
    eq_coefs = zeros(M,3);
    ro = .5;
    for i=1:M
        EQi = 1/ch_resp_est(i);
        EQ_min = 1/ch_resp_est(mod(i-1-1,M)+1);
        EQ_plu = 1/ch_resp_est(mod(i,M)+1);
        eqs = [EQ_min EQi EQ_plu];
        if eq_select == 2
            %   % geometric interpolation proposed by Aurelio & Bellanger
            %   % the coefficient computation from same paper
            %   % it's also approach 2 section 4.1.2 from D3.1
            EQ1 = EQ_min*(EQi/EQ_min)^ro;
            EQ2 = EQi*(EQ_plu/EQi)^ro;
         
        elseif eq_select == 3
            %   % approach 1 D3.1 section 4.1.1 linear interpolation
            EQ1 = (EQ_min+EQi)/2;
            EQ2 = (EQi+EQ_plu)/2;            
        end
        
        eq_coefs(i,1)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)+j*(EQ2-EQ1))/4;
        eq_coefs(i,2)= (EQ1+EQ2)/2;
        eq_coefs(i,3)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)-j*(EQ2-EQ1))/4; 
    
%         prim = upsample(conv(eq_coefs(i,:),scp_data(i,1:2:end-1)),2);
%         seco = circshift((upsample(conv(eq_coefs(i,:),scp_data(i,2:2:end)),2)).',1).';
%         psps = prim+seco;    
%         
    % %     re re ... im im ...
    %     sp_output(i,:) = psps(2:2*num_symbols+1);
    
    % %     re im re im ...
        equalizer_output(i,:)=conv(eq_coefs(i,:),scp_data(i,:));
    
        %sample
        sp_output(i,:) = equalizer_output(i,2*K+1+1:2*K+2*num_symbols+1);
    end
elseif eq_select == 4 %no equalizer
    for i=2*K+1:2*K+2*num_symbols
        sp_output(:,i-2*K) = scp_data(:,i);
    end
end

disp('+Subchannel processing is done.');
%% OQAM_Postprocessing
%

oqam_demod = zeros(M,2*num_symbols); %the matrix that will store demodulated signal
oqam_input = zeros(M,2*num_symbols); %this will hold after taking real part

mjn= (-j).^((1:2*num_symbols)-1); % (-j)^n

for k=1:M
    
    theta = mjn*(-j)^(k-1);
    oqam_input(k,:) =  real(sp_output(k,:).*theta);
    if ~mod((k-1),2) % even k
%         oqam_input(k,:) = real(sp_output(k,:).*theta);
        oqam_demod(k,:) = oqam_input(k,:).*repmat([1 j],1,num_symbols);
    else % odd k
%         oqam_input(k,:) = real(sp_output(k,:).*theta); 
        oqam_demod(k,:) = oqam_input(k,:).*repmat([j 1],1,num_symbols);
    end
end

disp('+OQAM Postprocessing is done.');
%% Symbol_Estimation
%

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
disp('+Symbol Estimation is done.');
%% Results

[error1,rate2] = symerr(m,m_est);
err=abs(m-m_est);
% dif_qam =[qam_m qam_est];
% error2 = symerr(qam_m(:,4:end),qam_est);
% dif_oqam =[oqam_m oqam_demod];
% error3 = symerr(oqam_m,oqam_demod);
% dif_m = [m m_est];
[error4,rate] = symerr(bits,bits_est);

disp(sprintf('Number of errorneous samples: %d/%d', error1,num_symbols*M))
disp(sprintf('SER: %f', rate2));
disp(sprintf('Number of bit errors: %d/%d', error4,num_bits))
disp(sprintf('BER: %f\n', rate))
disp('+Calculations and Statistics..');