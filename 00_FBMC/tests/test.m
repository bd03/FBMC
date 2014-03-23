%% test
%
% Burak Dayi
%
% This script will implement synthesis and analysis filterbanks
%
% Notes:
% Configuration will allow only subcarrier sizes of a power of two so as to
% be able to use FFT in the implementation
%
% The prototype filter used here is the PHYDYAS filter described in
% Deliverable 5.1: Prototype filter and structure optimization
%
% The implementation of the blocks follows what is described in Deliverable
% 5.1 transmultiplexer architecture.
%
% Channel implementation
% Subchannel processing
% Calculations and statistics will be handled later.
%
% Last updated: 06-03-2014

close all
clear all
clc
fprintf('--------------\n-----TEST-----\n--------------\n\n');

%% Configuration
%
% The user may set number of subcarriers, 
% overlapping factor, length of filter, prototype filter design parameters,
% modulation type, (number of frames, number of symbols), SNR value (and
% number of trials).
%
% The user input should be controlled and if
% necessary, error messages should be raised.

% -------------------------
% General filter parameters
% -------------------------
K = 4; % overlapping factor
M = 256; % number of subcarriers
% num_frames = 0; % number of frames
num_symbols = 100; % number of symbols
num_samples = M; %number of samples in a vector
modulation = 64; % 4-, 16- or 64-QAM
lp = K*M-1; % filter length
delay = K*M+1-lp; %delay requirement


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

if ~(ismember(modulation,[4 16 64 128 256]))
    modulation
    %error('Only 4-QAM scheme is supported. Please set modulation = 4');
    error('Only 4-,16-,64-QAM schemes are supported. Define the number of constellation points (4,16,64) with modulation variable.')
end

% --------------------------------
% Channel & simulation parameters
% --------------------------------
SNR = 1; %SNR value(s) in dB.
num_trials = 1; % number of trials desired

if num_trials<1
    num_trials
    error('Number of trials should be a positive integer.');
end

% ---------------------------------------
% Prototype filter frequency coefficients
% ---------------------------------------

% !This part should not be altered!

% K will select the row, the last column is background noise power in dB 
%(that might be come in handy in future) for reference
P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];

disp('+Configuration is obtained.');

%% Prototype_filter
%
% The following will compute the prototype filter coefficients and then plot 
% the prototype filter in frequency and time domain.
%
% Dependencies: lp, K, M, delay, P
% Output: h - prototype filter coefficients in time domain

disp('Prototype filter design')

h = zeros(1,lp); %coefficient array created with filter length
for m=1:lp
    h(m) = P(K,1);
    for k = 1:(K-1)
        h(m) = h(m) + 2*(((-1)^k)*P(K,k+1)*cos((2*pi*k*m)/(K*M)));
    end
end

%normalization of coefficients -> Do we really need this??
denom = P(K,1)+2*sum(P(K,2:4));
h=h/denom; %normalization

%delay implementation as on p.19 of deliverable 5.1
if delay>1
    h = [0 h 0];
    delay = delay - 1;
end

figure(11);
subplot(311);
plot(h);
title(sprintf('Time response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
%figure(12);
grid on;
subplot(312);
plot(abs(fft(h)));
title(sprintf('Freq response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
ylabel('Magnitude');
grid on;
subplot(313);
semilogy(abs(fft(h)));
title(sprintf('Freq magnitude response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
% ylabel('Log Magnitude');
grid on;
close all

disp('+Prototype filter is designed.');

%% Symbol_Creation
%
% This will create FBMC symbols with the modulation scheme defined
% in configuration file.
%
% Dependencies: modulation, num_subsymbol
% Output: qam_m - QAM modulated message

disp('Symbol Creation')

m = (randi(modulation,num_symbols,num_samples)-1).'; %random samples generated
qam_m = qammod(m, modulation); %qam modulation provided by MATLAB is used

disp('+Symbols are created.');

%% OQAM_Preprocessing
%
% This will perform OQAM preprocessing as it was described in
% PHYDYAS Deliverable 5.1
%
% Dependencies: qam_m - QAM modulated message 
% Output: oqam_m - OQAM modulated message

disp('OQAM preprocessing')

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

disp('+OQAM Preprocessing is done.');

%% Transmitter
%
% This will perform transmitter block (synthesis filter bank).
%
% Dependencies: oqam_m - OQAM modulated message
% Output: y - composite signal output

disp('Transmitter Block')
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
%transmitterda sorun olmamasi lazim

%% Channel
%
% This will perform channel impairments on the signal transmitted
% from the synthesis filter bank block.
%
% Dependencies: y - composite signal output
% Output: y_ch - composite signal gone through channel response

disp('Channel')

% channel impulse response
resp = [1]; % no impairments at the moment

%%%%%%%%%y_ch = conv(y,resp);
y_ch = y;

disp('+Channel effects are applied.');

%% Receiver
%
% This will perform receiver block (analysis filter bank).
%
% Dependencies: y_ch - composite signal gone through channel response
% Output: rx_output - receiver block output signal

disp('Receiver Block')

% The downsamplers following delay chain will form a Mx(K+1) matrix with
% channel output. Therefore, the delay chain and downsampler are jointly
% implemented as a vector reshape function.

% the matrix that reshaped input samples would be stored in
receiver_input = zeros(M,2*(K+num_symbols-1)); 

% reshaping (joint implementation of delay chain & downsamplers)
for r=0:K+num_symbols-1-1
%     y(1,1+r*M:M+r*M) = y(1,1+r*M:M+r*M)+ tx_poly_output(:,r+1).';
%     y(1,1+M/2+r*M:M+M/2+r*M) = y(1,1+M/2+r*M:M+M/2+r*M) + tx_poly_output(:,r+K+1).';
    receiver_input(:,r+1) = y(1,1+r*M:M+r*M);
    receiver_input(:,r+K+num_symbols) = y(1,1+M/2+r*M:M+M/2+r*M);    
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
    rx_fft_input(k,:) = [rx_poly_output(k,K:K+num_symbols-1) ...
        rx_poly_output(k, 2*K+num_symbols-2+K:2*K+num_symbols-2+K+num_symbols-1)];
    ppb(k,:) = b;
end

save('ppb.mat','ppb');

% fft performed
rx_fft_output=fft(rx_fft_input);

rx_output = zeros(M,2*num_symbols);

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
    
    for c = 1:num_symbols
        beta = [beta bb];
    end
   
    
    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    
    
    rx_output(k,1:2:2*num_symbols-1) =  rx_fft_output(k,1:num_symbols);
    rx_output(k,2:2:2*num_symbols) =  rx_fft_output(k,num_symbols+1:2*num_symbols);
    
    rx_output(k,:) = rx_output(k,:).*beta;
    
%     rx_output(k,:) = [];(1/sumfactor)*[sum(rx_fft_output(k,1:2*K-1)) ...
%         sum(rx_fft_output(k,2*K:2*(2*K-1)))].*beta; %(1/sumfactor)*
    
end

% rx_output will be sent to subchannel processing block

disp('+Receiver Block is processed.');

%% Subchannel_processing
%
% This script will perform subchannel processing steps.
%
% Dependencies: rx_output - receiver block output signal
% Output: sp_output - processed signal output

disp('Subchannel Processing')

sp_output = rx_output;

disp('+Subchannel processing is done.');

%% OQAM_Postprocessing
%
% This script will perform OQAM postprocessing as it was described in
% PHYDYAS Deliverable 5.1
%
% Dependencies: sp_output - processed signal output 
% Output: oqam_demod - estimated message

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

%% Symbol Estimation
%
% This will recombine the QAM samples and extract estimated symbol from
% those samples.
%
% Dependencies: oqam_demod - OQAM demodulated signal samples 
% Output: m_est - message estimation

qam_est = zeros(M,num_symbols);
m_est = zeros(M,num_symbols);
for k=1:M
    qam_est(k,:)= (1/sumfactor)*(oqam_demod(k,1:2:2*num_symbols-1)+oqam_demod(k,2:2:2*num_symbols));
    m_est(k,:) = qamdemod(qam_est(k,:),modulation);
end

disp('+Symbol Estimation is done.');
%% Results
error1 = symerr(m,m_est);
err=abs(m-m_est);
dif_qam =[qam_m qam_est];
error2 = symerr(qam_m,qam_est);
dif_oqam =[oqam_m oqam_demod];
error3 = symerr(oqam_m,oqam_demod);
dif_m = [m m_est];
disp(sprintf('Number of errors m/m_est: %d', error1))
disp(sprintf('Number of errors qam_dif: %d', error2))
disp(sprintf('Number of errors oqam_dif: %d', error3))