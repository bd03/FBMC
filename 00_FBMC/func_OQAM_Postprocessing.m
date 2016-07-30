function [oqam_demod_output] = func_OQAM_Postprocessing(sp_output,M,num_symbols,lower_boundary,upper_boundary)
%% func_OQAM_Postprocessing
%
% Burak Dayi
%
% This function will implement OQAM postprocessing
%
% Created: 24-11-2014

%% OQAM_Postprocessing

oqam_demod = zeros(M,2*num_symbols); %the matrix that will store demodulated signal
oqam_input = zeros(M,2*num_symbols); %this will hold after taking real part

mjn= (-j).^((1:2*num_symbols)-1); % (-j)^n

for k=lower_boundary:upper_boundary    
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

oqam_demod_output = oqam_demod(lower_boundary:upper_boundary,:);
% disp('+OQAM Postprocessing is done.');