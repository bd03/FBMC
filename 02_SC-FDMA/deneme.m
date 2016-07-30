%% deneme
clear all
close all
clc

M = 256;
N = 12;

num_sym=1000;
qam = 4;
num_bits = log2(qam)*num_sym*N;
mapping = 0;
Q=5;
start = 5;
normalization = [2 10 42 170];

bits = randi(2,1,num_bits)-1;

m = reshape(bi2de(reshape(bits,log2(qam),N*num_sym).','left-msb'),N,num_sym); % for now the array size is N*num_sym

qam_m = qammod(m, qam, pi/2,'gray')/sqrt(normalization(log2(qam)/2)); %built-in MATLAB qam qam 

tx_dft_output = fft(qam_m)/sqrt(N); % N-point DFT first.
tx_mapped = subcarrier_mapping(tx_dft_output,Q,M,1,0); % Subcarrier mapping
tx_ifft_output = ifft(tx_mapped)/sqrt(N/(M^2)); % M-point IFFT

rx_fft_output = fft(tx_ifft_output)/sqrt(M);%/sqrt(N/(M^2));%/sqrt(N);
rx_demapped = subcarrier_demapping(rx_fft_output,Q,N,1,0)/sqrt(M/N); % Subcarrier demapping
rx_output = ifft(rx_demapped)*sqrt(N);%/(N/M); % N-point DFT first.

m_est = zeros(N,num_sym);
for k=1:N
    m_est(k,:) = qamdemod(sqrt(normalization(log2(qam)/2))*rx_output(k,:),qam,pi/2,'gray');
%     m_est(k,:) = pskdemod(rx_output(k,:),qam,pi/2,'gray');
end

bits_est=reshape(de2bi(reshape(m_est,num_sym*N,1),'left-msb').',1,num_bits);

[error1,rate2] = symerr(m,m_est);
err=abs(m-m_est);
dif_m = [m m_est];
[error4,rate] = symerr(bits,bits_est);

disp(sprintf('Number of errorneous samples: %d/%d', error1,num_sym*M))
disp(sprintf('SER: %f', rate2));
%     disp(sprintf('Number of errors qam_dif: %d', error2))
%     disp(sprintf('Number of errors oqam_dif: %d', error3))
disp(sprintf('Number of bit errors: %d/%d', error4,num_bits))
disp(sprintf('BER: %f\n', rate))