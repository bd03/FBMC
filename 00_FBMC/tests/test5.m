% test5.m
%% implementation of system on p.13 of the journal with deconv

close all
%clear all
%clc

load('ppa.mat')
load('ppb.mat')
s= [1 -1;-1 1 ; 1 1; 1 1];

ifout =ifft(s);

%sfb ppt filter
for i=1:4
    aout(i,:) =conv(ppa(i,:),ifout(i,:));
end

%overlap & p/s
y = zeros(1,4*4+(2*size(s,2)-1)*4/2);
for i=1:4+size(s,2)-1
    y(1,1+4*(i-1):4*(i-1)+4)=y(1,1+4*(i-1):4*(i-1)+4)+aout(:,i).';
end


% s/p
for i=1:4+size(s,2)-1
    rin(:,i)=y(1,1+4*(i-1):4*(i-1)+4);
end

%afb ppt filter
for i=1:4
    bout(i,:)=conv(rin(i,:),ppb(i,:));
    coef(i,:) = conv(ppb(i,:),ppa(i,:));
end

%bout(2,:)=-1*bout(2,:);

%sampled=bout;
%sampled=bout(:,4:4+size(s,2)-1)./.6863;

 
fout = fft(bout);
fcoef = fft(coef);



figure(1);
for i=1:4
    inp=fout(i,2:length(fout(i,:)));
    s_est(i,:)=deconv(inp,fcoef(i,2:6));
%     subplot(2,2,i);
%     plot(real(fout(i,:)));
%     title(sprintf('sc=%d',i));
%     %axis([1 size(s,2) -2 2]);
end

% s_est=real(fout);
% 
% dif_s =[s, s_est];
% serr = abs(s-round(s_est));






























% %this works well
% %without overlapping

% close all
% clear all
% clc
% 
% load('ppa.mat')
% load('ppb.mat')
% s= [1 -1; 1 1; -1 1; 1 -1];
% 
% ifout =ifft(s);
% 
% 
% for i=1:4
%     filter(i,:)=conv(ppa(i,:),ppb(i,:));
%     out(i,:)=conv(ifout(i,:),filter(i,:));
% end
% 
% fout = fft(out);
% 
% figure(1);
% for i=1:4
%     subplot(2,2,i);
%     plot(real(fout(i,:)));
%     title(sprintf('sc=%d',i));
% end
% 
% 
% %stop
% s_est=real(fout(:,4:5));
% 
% dif_s =[s, s_est];

% this work with problem on sc=2
% w/ overlapping and sampling after fft

% 
% close all
% clear all
% clc
% 
% load('ppa.mat')
% load('ppb.mat')
% s= [1 -1; 1 1; -1 1; 1 -1];
% 
% ifout =ifft(s);
% 
% %sfb ppt filter
% for i=1:4
%     aout(i,:)=conv(ifout(i,:),ppa(i,:));
% end
% 
% %overlap & p/s
% y = zeros(1,12);
% for i=1:5
%     y(1,1+2*(i-1):2*(i-1)+4)=y(1,1+2*(i-1):2*(i-1)+4)+aout(:,i).';
% end
% 
% 
% % s/p
% for i=1:5
%     rin(:,i)=y(1,1+2*(i-1):2*(i-1)+4);
% end
% 
% %afb ppt filter
% for i=1:4
%     bout(i,:)=conv(rin(i,:),ppb(i,:));
% end
% 
% 
% 
%  
% fout = fft(bout);
% 
% figure(1);
% for i=1:4
%     subplot(2,2,i);
%     plot(real(fout(i,:)));
%     title(sprintf('sc=%d',i));
% end
% 
% s_est=real(fout(:,4:5));
% 
% dif_s =[s, s_est];