% first, load BER and CONF files

clear berr;

for i=1:length(conf.SNR_val)
	for k=1:length(conf.mod_sch)
		berr(i,k) = mean(BER(:,k,i));
	end
end

color_array={'y', 'm', 'c', 'r', 'g', 'b', 'k'};

figure;
semilogy(conf.SNR_val,berr(:,1), color_array{1+mod(1,length(color_array))})
hold on;

for i=2:length(conf.mod_sch)
	semilogy(conf.SNR_val,berr(:,i), color_array{1+mod(i,length(color_array))})
end

bertheory = berawgn(0:1:100,'qam',4);
semilogy(0:1:100,bertheory)
bertheory = berawgn(0:1:100,'qam',16);
semilogy(0:1:100,bertheory)
bertheory = berawgn(0:1:100,'qam',64);
semilogy(0:1:100,bertheory)
bertheory = berawgn(0:1:100,'qam',256);
semilogy(0:1:100,bertheory)
axis([0 SNR 1e-6 1])
grid on