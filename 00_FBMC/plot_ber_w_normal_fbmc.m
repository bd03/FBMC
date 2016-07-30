% first, load BER and CONF files

clear berr;

EbN0 = conf.sim_params.s_arr; % - round(10*log10(conf.channel.bw));
start_index = find(EbN0(1:end) == 0); % start from 0
if (size(start_index,2))==0
	start_index = 1;
end
% start_index = find(EbN0(1:end) == conf.sim_params.s_arr(1));

for i=1:length(conf.sim_params.s_arr)
	for k=1:length(conf.sim_params.q_arr)
		berr(i,k) = mean(BER(:,k,i));
	end
end

color_array={'y', 'm', 'c', 'r', 'g', 'b', 'k'};

figure;
% semilogy(0:1:(size(berr,1)-2),berr(start_index:end,1), color_array{1+mod(1,length(color_array))})
semilogy(EbN0(start_index:end),berr(start_index:end,1), color_array{1+mod(1,length(color_array))})
% semilogy(conf.sim_params.s_arr(start_index:end),berr(start_index:end,1), color_array{1+mod(1,length(color_array))})
hold on;

for i=2:length(conf.sim_params.q_arr)
	semilogy(EbN0(start_index:end),berr(start_index:end,i), color_array{1+mod(i,length(color_array))})
	% semilogy(conf.sim_params.s_arr(start_index:end),berr(start_index:end,i), color_array{1+mod(i,length(color_array))})
end

% figure;
% semilogy(conf.sim_params.s_arr,berr(:,1), color_array{1+mod(1,length(color_array))})
% hold on;

% for i=2:length(conf.sim_params.q_arr)
% 	semilogy(conf.sim_params.s_arr,berr(:,i), color_array{1+mod(i,length(color_array))})
% end

if (~conf.channel.fading)
	bertheory = berawgn(EbN0(1:end),'qam',4);
	semilogy(EbN0(start_index:end),bertheory(start_index:end),'--')
	bertheory = berawgn(EbN0(1:end),'qam',16);
	semilogy(EbN0(start_index:end),bertheory(start_index:end),'--')
	bertheory = berawgn(EbN0(1:end),'qam',64);
	semilogy(EbN0(start_index:end),bertheory(start_index:end),'--')
	bertheory = berawgn(EbN0(1:end),'qam',256);
	semilogy(EbN0(start_index:end),bertheory(start_index:end),'--')
end
axis([0 EbN0(end) 1e-6 1])
% axis([EbN0(1) EbN0(end) 1e-6 1])
grid on

xlabel('Normalized Eb/N0')
ylabel('BER')