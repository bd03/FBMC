                   
function [ci_imp_out] = LTE_channels (type,bandwidth)
% function [delay_a pow_a] = LTE_channels (type,bandwidth)
%LTE channels


% % EPA = 0;
% % ETU = 1;
% % EVA = 0;
% % 

bandw = bandwidth; % 5MHz

if type == 'EPA' % Low selectivity
     ci_imp = zeros(1,127);
     delay_a = [0 30 70 80 110 190 410]*1e-9;
     pow_a = [0 -1 -2 -3 -8 -17.2 -20.7];
elseif type == 'EVA' % Moderate selectivity
     ci_imp = zeros(1,127);
     delay_a = [0 30 150 310 370 710 1090 1730 2510 ]*1e-9;
     pow_a = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12-0 -16.9];
elseif type == 'ETU' % High selectivity
     ci_imp = zeros(1,127);
     delay_a = [0 50 120 200 230 500 1600 2300 5000]*1e-9;
     pow_a = [-1 -1.0 -1.0 0 0 0 -3 -5 -7];

else
    error('Invalid channel profile selection');
end


pow_a_lin = 10.^(pow_a./10);
% 
% %Making the sampled channel
tss = 1./bandw;
pos_a = round(delay_a./tss);
c_imp_sampled = [];
for i = min(pos_a):max(pos_a)
  c_imp_sampled(i+1) = sum(pow_a_lin(pos_a==i));
end
  ci_imp_out = sqrt((c_imp_sampled.^2./sum(c_imp_sampled.^2)));
