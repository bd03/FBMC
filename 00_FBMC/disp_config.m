if exist('mode_of_operation')
	if strcmp(mode_of_operation, 'SIMULATION')
		disp('Configuration:')
	    filter_params
	    supported_vals
	    sim_params
	    channel
	    sp_params
	elseif strcmp(mode_of_operation, 'MAIN')
		disp('Configuration:')
	    filter_params
	    sim_params
	    channel
	    sp_params
	end
elseif exist('conf')
	% that means it is a simulation result
	disp('Configuration:')
    conf.filter_params
    conf.supported_vals
    conf.sim_params
    conf.channel
    conf.sp_params
else
	error('No clue to print configuration.');
end