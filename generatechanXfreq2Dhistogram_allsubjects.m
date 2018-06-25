
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
relcnts_CF_subj = cell(length(InfantID),1);
BINrelcnts_CF_subj = cell(length(InfantID),1);
nTrlsPerSubj = cell(length(InfantID),1);

%% Extract artifact detected array for all subjects
subjcnt = 0;        cellarryemptynomoreflag = 0;
if ~isempty([relcnts_CF_subj{:,1}]), cellarryemptynomoreflag = 1; end

for ii=1:length(InfantID)
    %% Setting paths....
    disp(['Opening folder of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    serverPath1 = [InfantDir,infantfolder];
    cd(serverPath1)

    %% Load movement artifact identification array (w/trials)
        % Check if the EEG file exists
    fullFileName1 = 'importantWTCfreqs_TrialbyTrial_phaseextension_zero2pi.mat';   %correction.mat';   % important wavelet coherence events

    if ~exist(fullFileName1, 'file')
        % File does not exist.
        warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName1);
        disp(warningMessage)
        disp('Skipping to next infant data set')
        continue
    elseif cellarryemptynomoreflag == 1
        wtcfeatlist_summedalltrls = relcnts_CF_subj{ii,1};
        subjcnt = subjcnt + 1;
    else
        load(fullFileName1);
        subjcnt = subjcnt + 1;

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
    
        %% Save into cell array per subject
        nTrlsPerSubj{ii} = total_NoTrls;
        relcnts_CF_subj{ii,1} = wtcfeatlist_summedalltrls;
        relcnts_CF_subj{ii,2} = InfantID{ii};
        BINrelcnts_CF_subj{ii,1} = wtcfeatlist_summedalltrls > 0; % Binarized version
        BINrelcnts_CF_subj{ii,2} = InfantID{ii};    
    end
    %% Plot frequency-to-channel 2-D histogram per subject
    %
    figure(22)
    colormap(flipud(colormap('gray')))
    subtightplot(5,10,subjcnt,[.03 .01],[.1 .1],[.04 .01])
%     subplot(5,10,subjcnt)
    imagesc((1:64),log2(period_globalvar),wtcfeatlist_summedalltrls(:,GCN_idcs)*100);
    grid on, axis normal
    set(gca,'FontSize',14,'FontWeight','bold','FontName','Times New Roman',...
        'XMinorTick','on','YMinorTick','on')  
    caxis([0 20])    
    title([InfantID{ii},' (',num2str(nTrlsPerSubj{ii}),' trials)'],'fontsize',14,'fontweight','bold')

    if subjcnt == 45, xlabel('EEG Channel Regions','FontSize',18,...
    'FontWeight','bold','FontName','Times New Roman'), end
    if subjcnt == 21, ylabel('Frequency (Hz)','FontSize',18,...
    'FontWeight','bold','FontName','Times New Roman'), end

    if any(subjcnt == (41:49)),
        set(gca,'XTick',CGLs_idcs,'XtickLabel',(1:10))
%         rotateXLabels( gca(), 45 ),
    else
        set(gca,'XTick',CGLs_idcs,'XtickLabel',[]),
    end

    if any(subjcnt == (1:10:50)),
        yvalues = 1./period_globalvar;
        Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
        set(gca,'YTick',log2(1./Yticks(:)))
        set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
            'YDir','reverse', 'layer','top', ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
    else
        set(gca,'YTick',log2(1./Yticks(:)),'YtickLabel',[]),
    end
    %}   
end
cbh = colorbar('Location','north','Position',[0.919 0.109 0.055 0.019]);
ylabel(cbh,['% Artifactual',10,'Trials Detected'],'FontSize',14,...
    'FontWeight','bold','FontName','Times New Roman')
set(cbh,'ytick',(0:10:20),'yticklabel',(0:10:20),'FontSize',12,...
    'FontWeight','bold','FontName','Times New Roman')
myfig = figure(22);
set(myfig,'position',[1 1 1781 829]);

%% Save as PNG file
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'freqXchan2Dhistogram_allsubjects_v4.tif';
set(gcf,'PaperPositionMode','auto')
print(myfig,'-dpng', filename, '-r300')
