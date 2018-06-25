function [ V, COI, V_log, COI_log ] = plotCOIregion( N, M, fs )
%UNTITLED2 Summary of this function goes here
%   INPUTS:
%         N   --> number of time samples across both signals
%         M   --> number of frequency bins across both signals
%         fs  --> sampling frequeny (how many samples in one second)
%   OUTPUTS:
%         COI --> the cone of influence (COI) region shaded over a area
%                 which would make up the coherence plot.

Fourier_factor = 1.03; % based on Morlet wavelet Fourier wavelength (w_0 = 6)
e_foldtime = sqrt(2);  % based on Morlet wavelet e-folding time constant

coi_factor = Fourier_factor/e_foldtime;

t = (1:N);

% compute COI line
V = -((coi_factor/fs) * ( 1E-5 + ((N-1)/2) - abs(t - ((N-1)/2+1)) ));
V_log = -log((coi_factor/fs) * ( 1E-5 + ((N-1)/2) - abs(t - ((N-1)/2+1)) ));
freq = linspace(min(V),max(V),N/2);
freq_log = linspace(min(V_log),max(V_log),N/2);

% compute array where time-frequency elements are within the COI
COI=zeros(N/2,N);
COI_log=zeros(N/2,N);
for f = 1:N/2
    COI(f,:)= freq(f) <= V;
    COI_log(f,:)= freq_log(f) <= V_log;
end

%{
figure(1)
% plot COI line in a linear fashion
subplot(2,2,1)
plot(t,V,'r','linewidth',2)
title('COI Line Function')
ylabel('Frequency (Linear Scale)','fontweight','bold')
xlim([1 N])

subplot(2,2,2)
colormap(flipud(colormap('gray')))
imagesc(flipud(COI))  %fill the area that's inside the conditions
title('COI Region (shaded black)')
set(gca,'ytick',[],'yticklabel',[])

% plot COI line in a log fashion
subplot(2,2,3)
plot(t,V_log,'r','linewidth',2)
ylabel('Frequency (Log Scale)','fontweight','bold')
xlim([1 N]), ylim([min(V_log)-0.05 max(V_log)-10])

subplot(2,2,4)
colormap(flipud(colormap('gray')))
imagesc(flipud(COI_log))  %fill the area that's inside the conditions
xlim([1 N]), ylim([1+454 N/2+5])
set(gca,'ytick',[],'yticklabel',[])

% save
% myfig=figure(1);
% print(myfig,'-dpng','-r600','COIexample.png')
%}


figure(2)
% plot COI line in a log fashion
colormap(flipud(colormap('gray')))
imagesc(flipud(COI_log.*0.5))  %fill the area that's inside the conditions
xlim([1 N]), ylim([805 N/2]), caxis([0 1])
set(gca,'xtick',([1,500:500:2000]),'xticklabel',(-2:2))
set(gca,'ytick',(805:35:980),'yticklabel',fliplr(pow2(0:5)))
ylabel('Frequency (Hz)','fontweight','bold','fontsize',16,'fontname','times new roman')
xlabel('Time (s)','fontweight','bold','fontsize',16,'fontname','times new roman')
hold on
[~,contr_hndl] = contour(flipud(COI_log));
set(contr_hndl,'linewidth',2)
set(gca,'fontweight','bold','fontsize',14,'fontname','times new roman')
set(gcf,'position',[1 1 887 300])
save
myfig=figure(2);
print(myfig,'-dpng','-r600','COIexample_v2.png')

end

