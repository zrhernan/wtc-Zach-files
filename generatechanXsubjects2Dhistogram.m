
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

%% Preallocate
relcnts_CS_subj = nan(64,1);
CS_subjlabel = cell(1,1);
BINrelcnts_CS_subj = nan(64,1);
nTrlsPerSubj = cell(1,1);

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

        %% Sum across classes and frequency values
        total_NoTrls = sum(wtcfeatlist_trlSize);
        wtcfeatlist_summedtrls_mat = cat(3,wtcfeatlist_summedtrls{:});
        wtcfeatlist_summedalltrls = sum(sum(wtcfeatlist_summedtrls_mat,3),1)./(total_NoTrls*NoFreqBins);       

        %% Save into cell array per subject
        nTrlsPerSubj{subjcnt} = total_NoTrls;
        relcnts_CS_subj(:,subjcnt) = wtcfeatlist_summedalltrls;
        CS_subjlabel{subjcnt} = InfantID{ii};%[InfantID{ii},' [',num2str(total_NoTrls),' trials]'];
        BINrelcnts_CS_subj(:,subjcnt) = wtcfeatlist_summedalltrls > 0; % Binarized version
 
end

%% Display channel-to-subject 2-D histogram
myfig = figure(33);
subtightplot(1,4,(2:4),[.03 .01],[.15 .1],[.1 .05]);
colormap(flipud(colormap('gray')))
imagesc((1:length(CS_subjlabel)),(1:64),relcnts_CS_subj*100);
grid on, caxis([0 10]), axis normal
set(gca,'XTick',(1:length(CS_subjlabel)),'XtickLabel',CS_subjlabel)
set(gca,'YTick',CGLs_idcs,'YtickLabel',[],'YMinorTick','on')%(1:10)) %
set(gca,'FontSize',12,'FontWeight','bold','FontName','Times New Roman')

% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('Subject ID','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
% ylabel('EEG Channel Regions','FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
rotateXLabels( gca(), 45 )

cbh = colorbar;
ylabel(cbh,'% Trials Detected','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
set(cbh,'ytick',(0:2:10),'yticklabel',(0:2:10))
%{
%% Save as PNG file
set(myfig,'position',[1 1 1474 513]);
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'chanXsubject2Dhistogram.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r600')
%}
%%  Display channel 1-D histogram across all subjects
% myfig = figure(31);
subtightplot(1,4,1,[.03 .01],[.15 .1],[.1 .05]);
colormap(flipud(colormap('gray')))
bar((1:64),100*(sum(relcnts_CS_subj,2)/length(CS_subjlabel)),'FaceColor','k','BarWidth',1);
grid on, view([-90 -90])
set(gca,'XLim',[0.5 64.5],'XTick',CGLs_idcs,'XtickLabel',(1:10),...
    'FontSize',12,'FontWeight','bold','FontName','Times New Roman')
set(gca,'YLim',[0 5],'Ytick',(0:5),'YTickLabel',(0:5),'XMinorTick','on')
% rotateXLabels( gca(), 45 )
% title([InfantID{ii},' [',num2str(nTrls),' trials]'],'FontSize',14,...
%     'FontWeight','bold','FontName','Times New Roman')
xlabel('EEG Channel Regions','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
ylabel('% Trials Detected','FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')

%% Save as PNG file
% set(myfig,'position',[1 1 450 450]);
set(myfig,'position',[1 1 1868 938]);
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
% filename = 'chan1Dhistogram_groupanalysis.tif';
filename = 'chanXsubject2Dhistogram_AND_chan1Dhistogram_groupanalysis.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r300')