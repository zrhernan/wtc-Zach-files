
%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
disp_flag = 0;
[ InfantDataAnnotList, InfantID ] = defineInfantFolders( InfantDir, disp_flag );

%% Re-Arrange Channels by Region
load('BrainVision_1020_64ChannelOrder.mat')
[~, ~, ChanGroupLabels, CGLs_idcs, GCN_idcs ] = groupchansbyregion(channelOrder);

%% Initialize period (will convert to freq later) values
load('wtcperiodglobalvariable_v2.mat')
load('wtcbehaviorlistglobalvariable.mat')

%% Preallocate
relcnts_FS_subj = nan(length(period_globalvar),1);
FS_subjlabel = cell(1,1);
BINrelcnts_FS_subj = nan(length(period_globalvar),1);
nTrlsPerSubj = cell(1,1);
FAbysubj = cell(1,length(InfantID));

%% Extract artifact detected array for all subjects
subjcnt = 0;        cellarryemptynomoreflag = 0;
for ii=1:length(InfantID)
    %% Setting paths....
    disp(['Opening folder of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    serverPath1 = [InfantDir,infantfolder];
    cd(serverPath1)

    %% Load movement artifact identification array (w/trials)
        % Check if the EEG file exists
    fullFileName1 = 'importantWTCfreqs_TrialbyTrial_v2.mat';   % important wavelet coherence events

    if ~exist(fullFileName1, 'file')
        % File does not exist.
        warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName1);
        disp(warningMessage)
        disp('Skipping to next infant data set')
        continue
    else
        load(fullFileName1);
        subjcnt = subjcnt + 1;
    end
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

        %% Sum across classes and channels
        total_NoTrls = sum(wtcfeatlist_trlSize);
        wtcfeatlist_summedtrls_mat = cat(3,wtcfeatlist_summedtrls{:});
        wtcfeatlist_summedalltrls = sum(sum(wtcfeatlist_summedtrls_mat,3),2)./(total_NoTrls*64);       

        %% Save into cell array per subject
        nTrlsPerSubj{subjcnt} = total_NoTrls;
        relcnts_FS_subj(:,subjcnt) = wtcfeatlist_summedalltrls;
        FS_subjlabel{subjcnt} = InfantID{ii};%[InfantID{ii},' [',num2str(total_NoTrls),' trials]'];
        BINrelcnts_FS_subj(:,subjcnt) = wtcfeatlist_summedalltrls > 0; % Binarized version
        
        %% Load PSD information
        load('Kinematics\FTfreqanalysisPSDs.mat','trialavg')
        FAbysubj{subjcnt} = trialavg;
end

%% Grand average across subjects
cfg = [];
psd_subjavg = ft_freqgrandaverage(cfg,FAbysubj{:});

%% Display frequency-to-subject 2-D histogram
myfig = figure(33);
subtightplot(1,5,(3:5),[.03 .01],[.15 .1],[.1 .05]);
colormap(flipud(colormap('gray')))
imagesc((1:length(FS_subjlabel)),log2(period_globalvar),relcnts_FS_subj*100);
grid on, caxis([0 10]), axis normal
set(gca,'XTick',(1:length(FS_subjlabel)),'XtickLabel',FS_subjlabel,...
    'FontSize',12,'FontWeight','bold','FontName','Times New Roman',...
    'YMinorTick','on')

% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('Subject ID','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
% ylabel('Frequency (Hz)','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')

xvalues = 1./period_globalvar;
Xticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
set(gca,'YTick',log2(1./Xticks(:)))
set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
    'YDir','reverse', 'layer','top')
set(gca,'YTickLabel',[])%num2str(Xticks'))
rotateXLabels( gca(), 45 )

cbh = colorbar;
ylabel(cbh,'% Trials Detected','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
set(cbh,'ytick',(0:2:10),'yticklabel',(0:2:10))

%{
%% Save as PNG file
set(myfig,'position',[1 1 1474 513]);
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'freqXsubject2Dhistogram.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r300')
%}
%%  Display channel 1-D histogram across all subjects
% myfig = figure(31);
subtightplot(1,5,2,[.03 .01],[.15 .1],[.1 .05]);
colormap(flipud(colormap('gray')))
bar(log2(period_globalvar),100*(sum(relcnts_FS_subj,2)/length(FS_subjlabel)),'FaceColor','k','BarWidth',1);
grid on, view([-90 90])
xvalues = 1./period_globalvar;
Xticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
set(gca,'XTick',log2(1./Xticks(:)))
set(gca,'XLim',[-5.66 0.494],...
    'XDir','reverse', 'layer','top')
set(gca,'XTickLabel',[])%num2str(Xticks'))
set(gca,'YLim',[0 5],'Ytick',(0:5),'YTickLabel',(0:5))
set(gca,'FontSize',12,'FontWeight','bold','FontName','Times New Roman',...
    'XMinorTick','on','YMinorTick','on')
% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% xlabel('Frequency (Hz)','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
ylabel('% Trials Detected','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')

%% plot Grand-averaged PSD too
subtightplot(1,5,1,[.03 .01],[.15 .1],[.1 .05]);
xvalues_1 = psd_subjavg.freq;
mean_subjpsds = psd_subjavg.powspctrm;
plot(log2(xvalues_1),mean_subjpsds,'k','linewidth',2)
grid on, view([-90 90])
set(gca,'FontSize',12,'FontWeight','bold','FontName','Times New Roman')
set(gca,'XLim',[-0.494 5.66],'YLim',[0 0.1])
xlabel('Frequency (Hz)','FontSize',14,'FontWeight','bold','FontName','Times New Roman')
ylabel('Power (A.U.)','FontSize',14,'FontWeight','bold','FontName','Times New Roman')
set(gca,'XTick',log2(2.^(-1:5)),'XTickLabel',(2.^(-1:5)))


%% Save as PNG file
set(myfig,'position',[1 1 2068 938]);
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'freqXsubject2Dhistogram_AND_freq1Dhistogram_groupanalysis.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r300')