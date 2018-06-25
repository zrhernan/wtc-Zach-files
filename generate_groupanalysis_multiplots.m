
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
NoFreqBins = length(period_globalvar);

%% Preallocate (for frequency-subject 1D and 2D histograms)
relcnts_FS_subj = nan(length(period_globalvar),1);
FS_subjlabel = cell(1,1);
BINrelcnts_FS_subj = nan(length(period_globalvar),1);
FAbysubj = cell(1,1);

%% Preallocate (for channel-subject 1D and 2D histograms)
relcnts_CS_subj = nan(64,1);
CS_subjlabel = cell(1,1);
BINrelcnts_CS_subj = nan(64,1);

%% Extract artifact detected array for all subjects
subjcnt = 0;        nTrlsPerSubj = zeros(1,1);
if ~isempty(relcnts_CS_subj) && ~isempty(relcnts_FS_subj),
%     return
end
for ii=1:length(InfantID)
    %% Setting paths....
    disp(['Opening folder of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    serverPath1 = [InfantDir,infantfolder];
    cd(serverPath1)

    %% Load movement artifact identification array (w/trials)
        % Check if the EEG file exists
    fullFileName1 = 'importantWTCfreqs_TrialbyTrial_phaseextension_zero2pi.mat'; %v2.mat';   % important wavelet coherence events

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
        
        % save number of trials per subject
        total_NoTrls = sum(wtcfeatlist_trlSize);
        nTrlsPerSubj(subjcnt) = total_NoTrls; 
        
%% +++++++++++++++++++ CHANNEL-TO-SUBJECT HISTOGRAMS +++++++++++++++++++ %%
        % Sum across classes and frequency values
        wtcfeatlist_summedtrls_mat = cat(3,wtcfeatlist_summedtrls{:});
        wtcfeatlist_summedalltrls = sum(sum(wtcfeatlist_summedtrls_mat,3),1)./(total_NoTrls*NoFreqBins);       

        % Save into cell array per subject
        relcnts_CS_subj(:,subjcnt) = wtcfeatlist_summedalltrls;
        CS_subjlabel{subjcnt} = InfantID{ii};%[InfantID{ii},' [',num2str(total_NoTrls),' trials]'];
        BINrelcnts_CS_subj(:,subjcnt) = wtcfeatlist_summedalltrls > 0; % Binarized version

%% ++++++++++++++++++ FREQUENCY-TO-SUBJECT HISTOGRAMS ++++++++++++++++++ %% 
        % Sum across classes and channels
        wtcfeatlist_summedtrls_mat = cat(3,wtcfeatlist_summedtrls{:});
        wtcfeatlist_summedalltrls = sum(sum(wtcfeatlist_summedtrls_mat,3),2)./(total_NoTrls*64);       

        % Save into cell array per subject
        relcnts_FS_subj(:,subjcnt) = wtcfeatlist_summedalltrls;
        FS_subjlabel{subjcnt} = InfantID{ii};%[InfantID{ii},' [',num2str(total_NoTrls),' trials]'];
        BINrelcnts_FS_subj(:,subjcnt) = wtcfeatlist_summedalltrls > 0; % Binarized version
        
        %% Load PSD information (for frequency group analysis)
        load('Kinematics\FTfreqanalysisPSDs.mat','trialavg')
        FAbysubj{subjcnt} = trialavg;        
        
end

%% Grand average across subjects
cfg = [];
psd_subjavg = ft_freqgrandaverage(cfg,FAbysubj{:});

%% ++++++++++++++++++++ FREQUENCY-TO-SUBJECT PLOTS +++++++++++++++++++++ %%
%
%%%% Display frequency-to-subject 2-D histogram
myfig = figure(44);
gap = [.01 .01]; marg_h = [.13 .01]; marg_w = [.04 .01];
subtightplot(2,5,(3:5),gap,marg_h,marg_w);
colormap(flipud(colormap('gray')))
imagesc((1:length(FS_subjlabel)),log2(period_globalvar),relcnts_FS_subj*100);
grid on, caxis([0 20]), axis normal
set(gca,'XTick',(1:length(FS_subjlabel)),'XtickLabel',[])%FS_subjlabel)
set(gca,'FontSize',18,'FontWeight','bold','FontName','Times New Roman',...
    'YMinorTick','on')

% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% xlabel('Subject ID','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% ylabel('Frequency (Hz)','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')

xvalues = 1./period_globalvar;
Yticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
set(gca,'YTick',log2(1./Yticks(:)))
ylim(log2([period_globalvar(1),period_globalvar(end-2)]))
set(gca,'YDir','reverse', 'layer','top')
set(gca,'YTickLabel',[])%num2str(Yticks'))
rotateXLabels( gca(), 45 )

% cbh = colorbar;
% ylabel(cbh,'% Trials Detected','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% set(cbh,'ytick',(0:2:10),'yticklabel',(0:2:10))

%%%%  Display channel 1-D histogram across all subjects
subtightplot(2,5,2,gap,marg_h,marg_w);
colormap(flipud(colormap('gray')))
bar(log2(period_globalvar),100*(sum(relcnts_FS_subj,2)/length(FS_subjlabel)),'FaceColor','k','BarWidth',1);
grid on, view([-90 90])
xvalues = 1./period_globalvar;
Xticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
set(gca,'XTick',log2(1./Xticks(:)),'XMinorTick','on')
xlim(log2([period_globalvar(1),period_globalvar(end-2)]))
set(gca,'XDir','reverse', 'layer','top')
set(gca,'XTickLabel',[])%num2str(Xticks'))
set(gca,'YLim',[0 10],'Ytick',(0:10),'YTickLabel',[])%(0:5))
set(gca,'FontSize',18,'FontWeight','bold','FontName','Times New Roman',...
    'YMinorTick','on')
% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% xlabel('Frequency (Hz)','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
% ylabel('% Trials Detected','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')

%%%% plot Grand-averaged PSD too
subtightplot(2,5,1,gap,marg_h,marg_w);
xvalues_1 = psd_subjavg.freq;
mean_subjpsds = psd_subjavg.powspctrm;
plot(log2(xvalues_1),mean_subjpsds,'k','linewidth',2)
grid on, view([-90 90])
set(gca,'FontSize',18,'FontWeight','bold','FontName','Times New Roman')
set(gca,'XLim',[-0.329 5.671],'YLim',[0 0.1])
xlabel('Frequency (Hz)','FontSize',24,'FontWeight','bold','FontName','Times New Roman')
ylabel(['ACC Grand-Averaged*',10,'Power (A.U.)'],'FontSize',22,...
    'FontWeight','bold','FontName','Times New Roman')
text(5.126,0.095,['*',num2str(sum(nTrlsPerSubj)),' Trials Across ',...
    num2str(length(FS_subjlabel)),' Subjects'],...
    'fontweight','bold','fontname','times new roman','fontsize',14)
set(gca,'XTick',log2(2.^(-1:5)),'XTickLabel',(2.^(-1:5)),'XMinorTick','on')
set(gca,'YTick',(0:0.02:0.1),'YTickLabel',{'',(0.02:0.02:0.1)},'YMinorTick','on')
%}
%%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%%%

%% +++++++++++++++++++++ CHANNEL-TO-SUBJECT PLOTS ++++++++++++++++++++++ %%
%
%%%% Display channel-to-subject 2-D histogram
subtightplot(2,5,(8:10),gap,marg_h,marg_w);
colormap(flipud(colormap('gray')))
imagesc((1:length(CS_subjlabel)),(1:64),relcnts_CS_subj*100);
grid on, caxis([0 20]), axis normal
set(gca,'XTick',(1:length(CS_subjlabel)),'XtickLabel',CS_subjlabel)
set(gca,'YTick',CGLs_idcs,'YtickLabel',[],'YMinorTick','on')%(1:10)) %
set(gca,'FontSize',18,'FontWeight','bold','FontName','Times New Roman')

% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('Subject ID','FontSize',22,...
    'FontWeight','bold','FontName','Times New Roman')
% ylabel('EEG Channel Regions','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
rotateXLabels( gca(), 45 )

% generating colorbar for plot
cbh = colorbar('location','west','position',[0.154 0.148 0.017 0.245]);
ylabel(cbh,['% Artifactual',10,'Trials Detected'],'FontSize',18,...
    'FontWeight','bold','FontName','Times New Roman')
set(cbh,'ytick',(0:5:20),'yticklabel',(0:5:20),'FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')

%%%%  Display channel 1-D histogram across all subjects
subtightplot(2,5,7,gap,marg_h,marg_w);
colormap(flipud(colormap('gray')))
bar((1:64),100*(sum(relcnts_CS_subj,2)/length(CS_subjlabel)),'FaceColor','k','BarWidth',1);
grid on, view([-90 -90])
set(gca,'XLim',[0.5 64.5],'XTick',CGLs_idcs,'XtickLabel',(1:10),...
    'FontSize',18,'FontWeight','bold','FontName','Times New Roman','XMinorTick','on')
set(gca,'YLim',[0 10],'Ytick',(0:10),'YTickLabel',(0:10)) %[])%
% rotateXLabels( gca(), 45 )
% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('EEG Channel Regions','FontSize',24,...
    'FontWeight','bold','FontName','Times New Roman')
ylabel('% Artifactual Trials Detected','FontSize',24,...
    'FontWeight','bold','FontName','Times New Roman')
%}
%%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%%%
set(myfig,'position',[1 1 2170 970])%2439 1070]); %won't resize, use select arrow to open properties dialog box

%% Save as PNG file
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
% filename = 'chan1Dhistogram_groupanalysis.tif';
filename = 'groupanalysis_allplots_v4.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r600')
