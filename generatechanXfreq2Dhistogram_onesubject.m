
%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
disp_flag = 0;
[ InfantDataAnnotList, InfantID ] = defineInfantFolders( InfantDir, disp_flag );

%% Setting paths....
ii = input('Enter index from selected subject from ''InfantID'' list: ');
disp(['Opening folder of subject ',InfantID{ii}])
infantfolder = InfantDataAnnotList{ii};
serverPath1 = [InfantDir,infantfolder];
cd(serverPath1)

%% Load movement artifact identification array (w/trials)
load('importantWTCfreqs_TrialbyTrial_phaseextension_zero2pi.mat');   %correction.mat');

%% Sum across trials (but keep classes)
numCls = length(listofwtcfeats);
wtcfeatlist_summedtrls = cell(6,1);
wtcfeatlist_trlSize = zeros(6,1);
for cls = 1:numCls
    numTrls = size(listofwtcfeats{cls},3);
    if isempty(listofwtcfeats{cls}), numTrls = 0; end
    for trl = 1:numTrls
        wtcfeatlist_summedtrls{cls,1} = sum(listofwtcfeats{cls},3);
    end
    wtcfeatlist_trlSize(cls,1) = numTrls;
end

%% Sum across classes
total_NoTrls = sum(wtcfeatlist_trlSize);
wtcfeatlist_summedtrls_mat = cat(3,wtcfeatlist_summedtrls{:});
wtcfeatlist_summedalltrls = sum(wtcfeatlist_summedtrls_mat,3)./total_NoTrls;

%% Re-Arrange Channels by Region
load('BrainVision_1020_64ChannelOrder.mat')
[ChanGroupsWithLabels, ~, ChanGroupLabels, CGLs_idcs, GCN_idcs ] = groupchansbyregion(channelOrder);

%% Initialize period (will convert to freq later) values
load('wtcperiodglobalvariable_v2.mat')
load('wtcbehaviorlistglobalvariable.mat')

%% Display frequency-to-channel 2-D histogram
myfig = figure(11);
subtightplot(1,5,(2:5),[.03 .01],[.15 .1],[.1 .1]);
colormap(flipud(colormap('gray')))
imagesc( (1:64),log2(period_globalvar),wtcfeatlist_summedalltrls(:,GCN_idcs)*100);
grid on, caxis([0 20]), axis normal
set(gca,'XTick',CGLs_idcs,'XtickLabel',(1:10),'FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman',...
    'XMinorTick','on','YMinorTick','on')

% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('EEG Channel Regions','FontSize',16,...
    'FontWeight','bold','FontName','Times New Roman')
% ylabel('Frequency (Hz)','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')

yvalues = 1./period_globalvar;
Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
set(gca,'YTick',log2(1./Yticks(:)))
set(gca,'YDir','reverse', 'layer','top')
ylim(log2([period_globalvar(1),period_globalvar(end-2)]))
set(gca,'YTickLabel',[])%num2str(Yticks'))
% rotateXLabels( gca(), 45 )
hold on
cbh = colorbar('position',[0.925 0.2 0.01 0.65]);
ylabel(cbh,'% Artifactual Trials Detected','FontSize',16,...
    'FontWeight','bold','FontName','Times New Roman')
set(cbh,'ytick',(0:2:20),'yticklabel',(0:2:20))

%% Load PSD information
load('Kinematics\FTfreqanalysisPSDs.mat')

%% Display the Magnitude Head Acceleration PSD
subtightplot(1,5,1,[.03 .01],[.15 .1],[.1 .1]);
xvalues = pertrial.freq;
psd_trls = squeeze(pertrial.powspctrm);
mean_psds = mean(psd_trls);    stnderr_psds = std(psd_trls,[],1);
LPs.col = {'k'};
mseb(log2(xvalues),mean_psds,stnderr_psds,LPs); view([-90 90]), grid on
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times New Roman')
xlim(-fliplr(log2([period_globalvar(1),period_globalvar(end-2)]))), ylim([0 0.1])
xlabel('Frequency (Hz)','FontSize',16,'FontWeight','bold','FontName','Times New Roman')
ylabel('ACC Power (A.U.)','FontSize',16,'FontWeight','bold','FontName','Times New Roman')
set(gca,'XTick',log2(2.^(-1:5)),'XTickLabel',(2.^(-1:5)))

    
%% Save as PNG file
set(myfig,'position',[1 1 752 581]);
    cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'freqXchan2Dhistogram_RB23_v4.tif';
% cd([serverPath1,'\Kinematics'])
print(myfig,'-dpng', filename, '-r300')

