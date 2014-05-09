%% zz_OneFileFBMC5.m
%
% Burak Dayi

% Sending multiple frames with preambles of three oqam symbols.

% The collection of three estimation methods implemented in
% OneFileFBMC2-3 & 4

%
% Created: 08-05-2014

%close all
clear all
clc
fprintf('--------------\n-----FBMC-----\n--------------\n\n');

%% Transmission 
%% Config

% 2.b Main mode parameters
%---- General filterbank parameters ----%
K = 4; % overlapping factor 
M = 256; % number of subcarriers
num_frames = 10; % number of data frames in each FBMC block
syms_per_frame = 10; %number of symbols per FBMC frame
num_symbols = num_frames*syms_per_frame; % total number of data symbols
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

%---- Channel settings ----%
% noise settings
noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
SNR = 10; % SNR of the channel. noisy=1 to see the effects on channel

%rayleigh channel settings
fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
bw = 5e+6; % Transmission Bandwidth
fd= 0; %max. doppler shift
channel_profiles = ['EPA', 'EVA', 'ETU']; % Valid channel profile selections
profile ='ETU'; %Channel profile
[delay_a pow_a] = LTE_channels (profile,bw);
ch_resp = rayleighchan(1/bw,fd,delay_a,pow_a); %channel model
%ch_resp.storeHistory = 1;
ch_resp.storePathGains =1;

%---- Channel estimation settings ----%
% method
valid_methods = ['IAM', 'IAM4']; % 'POP' may be added later.
estimation_method = 'IAM';

% preamble
% IAM preambles 
preamble = [zeros(M,1) repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];
% preamble = [zeros(M,1) repmat([1 -1 -1 1].',M/4,1) zeros(M,1)];
% preamble = [zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1)];
% POP preambles
% preamble = [repmat([1 -1].',M/2,1) zeros(M,1)];
% IAM4 preambles
% preamble = [zeros(M,1) zeros(M,1) repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];
if strcmp(estimation_method,'IAM4')
    preamble =[zeros(M,1) preamble];
end


%---- Equalizer settings ----%

eq_select = 2; % selection of equalizer type 1: one tap, 
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

if ~ismember(valid_methods,estimation_method)
    estimation_method
    error('estimation_method should be either IAM, IAM4 or POP.');
end


% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if noisy
    disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, fd=%d, SNR=%d dB', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,fd,SNR));
else
    disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, fd=%d, SNR=Ideal', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,fd));
end
disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
%pause;

disp('+Configuration is obtained.');
Prototype_filter;
disp('+Prototype filter is designed.');
%% Symbol_Creation

% bit sequence creation
bits = randi(2,1,num_bits)-1;
% load('bits.mat');
%save('bits.mat','bits');

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

jp = zeros(size(oqam_m));

% preamble insertion
for i=num_frames:-1:1
    oqam_m = [oqam_m(:,1:2*(i-1)*syms_per_frame) preamble oqam_m(:,2*(i-1)*syms_per_frame+1:end)];
    jp = [jp(:,1:2*(i-1)*syms_per_frame) preamble jp(:,2*(i-1)*syms_per_frame+1:end)];
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
disp('+OQAM Preprocessing is done.');
%% Transmitter

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
ppa = zeros(M,K);
betas = zeros(M,num_oqam_subsymbols);

for k=1:M
    
    beta = ((-1).^((k-1)*(1:num_oqam_subsymbols)))*((-1).^((k-1)*K));
    
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
disp('+Transmitter Block is processed.');

%% Channel

if fading
    y_response = filter(ch_resp,ones(size(y)));
    y_filtered = filter(ch_resp,y);
    y_filtered2 = y_response.*y;
%     y_filtered = y_filtered2;
    y_filtered3 = conv(y_response,y);
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
receiver_input_1 = zeros(M,K+ceil(num_oqam_subsymbols/2)-1); 
receiver_input_2 = zeros(M,K+floor(num_oqam_subsymbols/2)-1);

% reshaping will be separated
% reshaping (joint implementation of delay chain & downsamplers)
% PPN1
for r=0:K+ceil(num_oqam_subsymbols/2)-1-1
    receiver_input_1(:,r+1) = y_ch(1,1+r*M:M+r*M);
end

for r=0:K+floor(num_oqam_subsymbols/2)-1-1
    receiver_input_2(:,r+1) = y_ch(1,1+M/2+r*M:M+M/2+r*M);
end

rx_poly_output_1 = zeros(M,K+ceil(num_oqam_subsymbols/2)-1+K-1); %output of PPN1 will be stored in this
rx_poly_output_2 = zeros(M,K+floor(num_oqam_subsymbols/2)-1+K-1); %output of PPN2 will be stored in this

jp_rx_poly_output = zeros(M,2*(K+num_symbols-1+K-1)); %output of polyphase filters will be stored in this
jp_rx_fft_input = zeros(M,2*num_symbols);

% polyphase filter coefficients for archive
ppb = zeros(M,K);

% polyphase filters are applied
for k=1:M
    b = h(M-k+1:M:lp+1); % related polyphase filter coefficients sieved
    % we use lp+1 b/c of the delay implemented in prototype filter design   
    rx_poly_output_1(k,:) = conv(receiver_input_1(k,:),b);
    rx_poly_output_2(k,:) = conv(receiver_input_2(k,:),b);
    ppb(k,:) = b;
end

% rearrangement
len=K+ceil(num_oqam_subsymbols/2)-1+K-1+K+floor(num_oqam_subsymbols/2)-1+K-1;
rx_fft_input = zeros(M,len);
rx_fft_input(:,1:2:end)=rx_poly_output_1;
rx_fft_input(:,2:2:end)=rx_poly_output_2;

save('ppb.mat','ppb');

% fft performed
rx_fft_output=fft(rx_fft_input);
jp_rx_fft_output=fft(jp_rx_poly_output);

rx_output = zeros(M,len);
jp_rx_output = zeros(M,2*(num_symbols+K-1+K-1));

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
    beta = ((-1).^((k-1)*(1:len)))*((-1).^((k-1)*K));

    % ifft input is oqam_m on kth subchannel multiplied by beta multiplier
    % on that subchannel, in accordance with the order of the OQAM
    % subsymbol
    
    rx_output(k,:) = rx_fft_output(k,:).*beta;
end

% rx_output will be sent to subchannel processing block
disp('+Receiver Block is processed.');

%% Subchannel_processing

scp_input = rx_output(:,2*K-1:end-(2*K-1)+1);

if strcmp(estimation_method,'IAM')
    % remove preamble
    for i=1:num_frames
        scp_preamble(:,i) = scp_input(:,2+(i-1)*(syms_per_frame+1)*2+(i-1));
        scp_data(:,1+(i-1)*(syms_per_frame)*2:i*(syms_per_frame)*2) = scp_input(:,4+(i-1)*(syms_per_frame+1)*2+(i-1):4+(i-1)*(syms_per_frame+1)*2+(i-1)+syms_per_frame*2-1);
    end
    
    % estimation of channel
    for i=1:num_frames
        ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(:,2)*sumfactor);
    end

    if eq_select == 1
        % one tap equalizer
        eq_coefs=[];
        for ii=1:num_frames
            eq_coefs =[eq_coefs repmat(ch_resp_est(:,ii),1,syms_per_frame*2)];
        end

        sp_output = scp_data./eq_coefs;
    elseif eq_select == 2 || eq_select == 3
        % % three tap equalizer
        % % coef computation
        eq_coefs = zeros(M,3);
        ro = .5;
        for i=1:M
            for ii=1:num_frames
                EQi = 1/ch_resp_est(i,ii);
                EQ_min = 1/ch_resp_est(mod(i-1-1,M)+1,ii);
                EQ_plu = 1/ch_resp_est(mod(i,M)+1,ii);
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

                % %     re im re im ...
                equalizer_output =conv(eq_coefs(i,:),rx_output(i,:));
                sp_output(i,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
                    = equalizer_output(2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1):2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1)+2*syms_per_frame-1);
            end
        end
    elseif eq_select == 4 %no equalizer
        sp_output=scp_data;
    end
elseif strcmp(estimation_method,'IAM4')
    % remove preamble
    for i=1:num_frames
        scp_preamble(:,i) = scp_input(:,3+(i-1)*(syms_per_frame+2)*2);
        scp_data(:,1+(i-1)*(syms_per_frame)*2:i*(syms_per_frame)*2) = scp_input(:,5+(i-1)*(syms_per_frame+2)*2:5+(i-1)*(syms_per_frame+2)*2+syms_per_frame*2-1);
    end
    
    %
    % %estimation of channel
    for i=1:num_frames
        ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(:,3)*sumfactor);
    end

    if eq_select == 1
        % one tap equalizer
        eq_coefs=[];
        for ii=1:num_frames
            eq_coefs =[eq_coefs repmat(ch_resp_est(:,ii),1,syms_per_frame*2)];
        end

        sp_output = scp_data./eq_coefs;
    elseif eq_select == 2 || eq_select == 3
        % % three tap equalizer
        % % coef computation
        eq_coefs = zeros(M,3);
        ro = .5;
        for i=1:M
            for ii=1:num_frames
                EQi = 1/ch_resp_est(i,ii);
                EQ_min = 1/ch_resp_est(mod(i-1-1,M)+1,ii);
                EQ_plu = 1/ch_resp_est(mod(i,M)+1,ii);
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

                % %     re im re im ...
                equalizer_output =conv(eq_coefs(i,:),rx_output(i,:));
                sp_output(i,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
                    = equalizer_output(2*K+2+1+(ii-1)*2*(syms_per_frame+2)+1:1+2+2*K+2*syms_per_frame+(ii-1)*2*(syms_per_frame+2));
            end
        end
    elseif eq_select == 4 %no equalizer
        sp_output=scp_data;
    end
end

disp('+Subchannel processing is done.');
%% OQAM_Postprocessing

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
disp('+Symbol Estimation is done.');

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
bit_err = abs(bits-bits_est);
% dif_qam =[qam_m qam_est];
% error2 = symerr(qam_m(:,4:end),qam_est);
% dif_oqam =[oqam_m oqam_demod];
% error3 = symerr(oqam_m,oqam_demod);
% dif_m = [m m_est];
[error4,rate] = symerr(bits,bits_est);
disp(sprintf('Number of errorneous samples: %d/%d', error1,num_symbols*M))
disp(sprintf('SER: %f', rate2));
%     disp(sprintf('Number of errors qam_dif: %d', error2))
%     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
disp(sprintf('Number of bit errors: %d/%d', error4,num_bits))
disp(sprintf('BER: %f\n', rate))

max_num_err=0;
for i=1:num_symbols
    num_errored(i) = length(find(err(:,i)));
    num_errored_bits_per_symbol(i) = length(find(bit_err(1+(i-1)*num_bits/num_symbols:i*num_bits/num_symbols)));
    
    if max_num_err<num_errored(i)
        max_num_err=num_errored(i);
    end
end    

errors = zeros(num_symbols,max_num_err);

for i=1:num_symbols
    errors(i,1:length(find(err(:,i))))= find(err(:,i));
end

for i=1:num_frames
    num_errored_bits_per_frame(i) = length(find(bit_err(1+(i-1)*num_bits/num_frames:i*num_bits/num_frames)));
end

figfig=figure(54);
subplot(311)
plot(1:num_symbols,num_errored)
xlabel('FBMC symbols');
ylabel('Number of errors')
title('Number of errors w.r.t symbols')
subplot(312)
plot(1:num_symbols,num_errored_bits_per_symbol)
xlabel('FBMC symbols');
ylabel('Number of errors')
title('Number of bit errors w.r.t symbols')
subplot(313)
plot(1:num_frames,num_errored_bits_per_frame)
xlabel('FBMC frames');
ylabel('Number of errors')
title('Number of bit errors w.r.t frames')

figu=figure(56);
plot(1:num_symbols,errors,'o')
xlabel('Symbols')
ylabel('Subcarrier')
title('Errors on each SC w.r.t symbols')
c1=clock;
% savefig(figu,sprintf('%d%d%d%d%2d%d-figu',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
% savefig(figfig,sprintf('%d%d%d%d%2d%d-figfig',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
% save(sprintf('%d%d%d%d%2d%d-ch_resp.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'ch_resp');
% save(sprintf('%d%d%d%d%2d%d-bits.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'bits');