%% func_preamble_creation: function description
function [preamble,length_preamble,est_col] = func_preamble_creation(M, preamble_sel, zero_pads, extra_zero, user_indices, eq_select, fractional)
%% func_Analysis_Filter_Bank
%
% Burak Dayi
%
% This function will return the preamble. 
%
% Created: 25-02-2015

preamble = NaN;
center_preamble = NaN;
est_col = 0;
length_preamble = 0;

% define preambles here
switch preamble_sel
	case 0
		center_preamble = repmat([1 -j -1 j].',M/4,1); % IAM-I
		est_col = 1+zero_pads; % estimation on this column
	case 1
		center_preamble = repmat([1 1 -1 -1].',M/4,1); % IAM-R
		est_col = 1+zero_pads; % estimation on this column
	case 2
		center_preamble = repmat(repmat([1 -j -1 j].',M/4,1),1,3); % IAM-I with triple repetition.
		est_col = 2+zero_pads; % estimation on middle column
	otherwise
		error('Invalid preamble selection.')
end

if fractional
	% now according to equalizer and user indices, 
	% the preamble will be prepared
	% the unused subchannels will be sieved out.
	if eq_select == 1 % 1: one tap
		% we need only those indices over which data is transmitted.
		if preamble_sel == 2
			center_preamble(1:(user_indices(1)-1),:) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+1):(user_indices(2*i+1)-1),:) = 0;
			end

			center_preamble((user_indices(end)+1):end,:) = 0;		
		else
			center_preamble(1:(user_indices(1)-1)) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+1):(user_indices(2*i+1)-1)) = 0;
			end

			center_preamble((user_indices(end)+1):end) = 0;
		end
	elseif (eq_select == 2) | (eq_select == 3) %2: % three taps
		% when three taps are applied, in order to get a reliable channel estimation
		% on edge frequencies, we need to extend the covered subchannels by one
		% from each upper and lower bound 

		if preamble_sel == 2
			center_preamble(1:(user_indices(1)-2),:) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+2):(user_indices(2*i+1)-2),:) = 0;
			end

			center_preamble((user_indices(end)+2):end,:) = 0;
		else	
			center_preamble(1:(user_indices(1)-2)) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+2):(user_indices(2*i+1)-2)) = 0;
			end

			center_preamble((user_indices(end)+2):end) = 0;
		end
	elseif eq_select == 4 % no equalizer
		% we need only those indices over which data is transmitted. ?!?!?!?!?!
		if preamble_sel == 2
			center_preamble(1:(user_indices(1)-1),:) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+1):(user_indices(2*i+1)-1),:) = 0;
			end

			center_preamble((user_indices(end)+1):end,:) = 0;		
		else
			center_preamble(1:(user_indices(1)-1)) = 0;

			for i=1:((length(user_indices)/2)-1)
				center_preamble((user_indices(2*i)+1):(user_indices(2*i+1)-1)) = 0;
			end

			center_preamble((user_indices(end)+1):end) = 0;
		end
	else
		error('Unhandled equalizer selection.')
	end
end

% construct preamble with zero pads
pads = repmat(zeros(M,1),1,zero_pads);
early_preamble = [pads center_preamble pads];
if extra_zero
	preamble = [early_preamble zeros(M,1)];
else
	preamble = early_preamble;	
end

length_preamble = size(preamble,1)*size(preamble,2);