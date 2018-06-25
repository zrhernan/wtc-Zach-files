%% Load file
%{
[~, name] = system('hostname');
if ismember(cellstr(name),{'Beezi11-SuperHP','BMI-Gait-2'})
    usr = 'zrhernan';
elseif ismember(cellstr(name),{'Zach-HPenvy17'})
    usr = 'Zachery';
end

load(['C:\Users\',usr,'\OneDrive\WTC_resultsbyTrial.mat'])
%}

%% Generate List of Infant Data Folders
InfantDir = '\\172.27.216.40\Contreras-UH\Infantdata\Infantdata\Data\';
disp_flag = 0;
[ InfantDataAnnotList, InfantID ] = defineInfantFolders( InfantDir, disp_flag );

%%
list_wtcfeats_Xsubjs = cell(length(InfantID),1);
numTrl_wtcfeats_Xsubjs = cell(length(InfantID),1);
tic;
for ii=1:length(InfantID)
    disp(['Opening folder of subject ',InfantID{ii}])
    infantfolder = InfantDataAnnotList{ii};
    serverPath1 = [InfantDir,infantfolder];
    cd(serverPath1)

    %% Load movement artifact identification array (w/trials)
    % Check if the EEG file exists
    fullFileName1 = 'importantWTCfreqs_TrialbyTrial_phaseextension_zero2pi.mat';   % important wavelet coherence events

    if ~exist(fullFileName1, 'file')
        % File does not exist.
        warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName1);
        disp(warningMessage)
        disp('Skipping to next infant data set')
        continue
    else
        load(fullFileName1);
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
    end
%{
    %% Generate binary list of frequencies match all criteria:
    % - are statistically significant (based on Monte Carlo tests)
    % - have squared magnitude coherence > 0.5
    % - is within the Cone of Influence
    % - phase is between 0 and -90 degrees (i.e. in-phase to EEG leading acceleration), and
    % - covers 10% of the total time range (dependent on frequency due to the Cone of Influence
    [ listofwtcfeats ] = detect_importantWTCs( WTCresults.Rsq, WTCresults.period, WTCresults.coi, WTCresults.wtcsig, WTCresults.Wxy );
    %
    listofwtcfeats = cell(size(WTCresults.Rsq,1),1);
    for cls = 1:size(WTCresults.Rsq,1)
        for chn = 1:size(WTCresults.Rsq,2)
            % Region of Contours where coherence values are statistically significant
            signifRegions = double(WTCresults.wtcsig{cls,chn} >= 1);

            % Region of Cone of Influence
            outsidecoi = nan(size(WTCresults.wtcsig{cls,chn}));
            for s=1:length(WTCresults.period{cls,chn})
                outsidecoi(s,:)=(WTCresults.period{cls,chn}(s)<=WTCresults.coi{cls,chn});
            end

            % Region of Phases where Phase is between 0 and -90 (i.e. where 
            aWxy=angle(WTCresults.Wxy{cls,chn});
            aaa=aWxy;
            aaa(WTCresults.Rsq{cls,chn} < 0.5) = NaN; %remove phase indication where Rsq is low
            leadingphaseRegions = aaa >= 0 & aaa < pi/2;

            % Find the Intersection of all 3 regions
            combineRegions = signifRegions & outsidecoi & leadingphaseRegions;

            % Indicate frequencies where the region covers 10% of the time or more 
            signif_10prcntoftime = zeros(size(combineRegions,1),1);
            time_thr = 10;
            for prd = 1:size(combineRegions,1)
                if sum(combineRegions(prd,:)) >= round(sum(outsidecoi(1,:))/time_thr)
                    signif_10prcntoftime(prd) = 1;
                end
            end

            listofwtcfeats{cls}(:,chn) = signif_10prcntoftime;
        end
    end
    %}
    list_wtcfeats_Xsubjs{ii} = wtcfeatlist_summedtrls;
    numTrl_wtcfeats_Xsubjs{ii} = wtcfeatlist_trlSize;
    
end
disp('Done finding important coherence values')
timecompute(toc);

%% Split by age group
consolid_list = squeeze(cat(3,list_wtcfeats_Xsubjs{:}));
agegroups = {{6},{7,8,9},{10,11,12},{13,14,15,16,17,18},{19,20,21,22,23,24}};
consolid_list_agegroup = cell(length(agegroups),1);
numTrl_list_agegroup = cell(length(agegroups),1);
subjlist_agegroup = cell(length(agegroups),1);
agelist = zeros(length(InfantID),1);
for ii=1:length(InfantID)
    if isempty(list_wtcfeats_Xsubjs{ii}), continue, end
    agelist(ii) = str2double(cellstr(InfantID{ii}(regexp(InfantID{ii},'\d'))));
    for ag = 1:length(agegroups)
        if any(agelist(ii) == [agegroups{ag}{:}])
            consolid_list_agegroup{ag} = [consolid_list_agegroup{ag},list_wtcfeats_Xsubjs{ii}];
            numTrl_list_agegroup{ag} = [numTrl_list_agegroup{ag},numTrl_wtcfeats_Xsubjs{ii}];
            subjlist_agegroup{ag}{end+1,1} = InfantID{ii};
        end
    end
end
%% Sum relevant coherence moments by age group
WTC_cnts_agegroup = cell(size(consolid_list_agegroup{1},1),length(agegroups));
WTC_max_agegroup = zeros(size(consolid_list_agegroup{1},1),length(agegroups));
WTC_numTrls_agegroup = zeros(size(consolid_list_agegroup{1},1),length(agegroups));
WTC_numSubjs_agegroup = zeros(length(agegroups),1);
for ag = 1:length(agegroups)
    for bhvr = 1:size(consolid_list_agegroup{1},1)
        WTC_cnts_agegroup{bhvr,ag} = sum(cat(3,consolid_list_agegroup{ag}{bhvr,~cellfun(@isempty,{consolid_list_agegroup{ag}{bhvr,:}})}),3);
        WTC_max_agegroup(bhvr,ag) = max(WTC_cnts_agegroup{bhvr,ag}(:));
        WTC_numTrls_agegroup(:,ag) = sum(numTrl_list_agegroup{ag},2);
    end
    WTC_numSubjs_agegroup(ag) = size(consolid_list_agegroup{ag},2);
end

%% Re-Arrange Channels by Region
load('BrainVision_1020_64ChannelOrder.mat')
[~, ~, ChanGroupLabels, CGLs_idcs, GCN_idcs ] = groupchansbyregion(channelOrder);

%{
%% Plot them
myfig = figure(11);
agegroup_str = {'6 mo.','7-9 mo.','10-12 mo.','13-18 mo.','19-24 mo.'};
load('wtcperiodglobalvariable_v2.mat')
load('wtcbehaviorlistglobalvariable.mat')
array2disp = WTC_cnts_agegroup;
cmax = max(WTC_max_agegroup(:));
cmin = 0;
numSProws = length(agegroups)+1;
numSPcols = size(array2disp,1)*7+2;

for ag = 1:length(agegroups)+1
    SProwIDX = 0;
    for cls = 1:size(array2disp,1)
  
        % check if data structure exists
        if isempty(array2disp{cls,ag}), continue, end

        %% plotting the heatmap
        rowshift = (ag-1)*(size(array2disp,1)*7+2);
        SProwIDX = SProwIDX + 1;        
        
        HMplot_beg = 7*(SProwIDX-1)+1+rowshift;       HMplot_end = HMplot_beg+4;
        subplot(numSProws,numSPcols,HMplot_beg:HMplot_end)
        colormap(flipud(colormap('gray')))
        imagesc( (1:64),log2(period_globalvar),array2disp{cls,ag}(:,GCN_idcs)./(WTC_numTrls_agegroup(cls,ag)) );%WTC_max_agegroup(cls,ag));
        grid on, caxis([cmin 0.2])
        set(gca,'XTick',CGLs_idcs,'XtickLabel',ChanGroupLabels,...
            'FontSize',12,'FontWeight','bold')

        yvalues = 1./period_globalvar;
        Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
        set(gca,'YTick',log2(1./Yticks(:)))
        set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
            'YDir','reverse', 'layer','top', ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
        rotateXLabels( gca(), 45 )
        
        % Add lines to emphasize frequencies of importance
        chIDX = (1:length(channelOrder));
        % for line at 1 Hz
        yval1 = repmat(log2(1),1,length(chIDX));  
        line(chIDX,yval1,'color','r','linestyle','--','linewidth',2)
        % for line at 12 Hz
        yval2 = repmat(log2(1/12),1,length(chIDX));  
        line(chIDX,yval2,'color','r','linestyle','--','Linewidth',2)     
        if ag==1, title(upper(ClassOrderList{cls}),'FontSize',14,'FontWeight','bold'), end
%         if any(spIDX == [27 28]), xlabel('Channel Groups'), end
        if cls==1, ylabel([agegroup_str{ag},10,10],'FontSize',14,'FontWeight','bold'),end
        if ag==3 && cls==1, ylabel([agegroup_str{ag},10,10,'Frequency (Hz)'],'FontSize',14,'FontWeight','bold'), end

        %% plotting the histogram across channels
        Hplot_idx = HMplot_end+1;
        subplot(numSProws,numSPcols,Hplot_idx)
        hist_chns = sum(array2disp{cls,ag},2)./(length(channelOrder)*WTC_numTrls_agegroup(cls,ag));
        barh(log2(period_globalvar),hist_chns,'FaceColor','k')
        yvalues = 1./period_globalvar;
        Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
        set(gca,'YTick',[])
        set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
            'YDir','reverse', 'layer','top')
        set(gca,'XLim',[0 0.2],'Xtick',(0:0.2:0.2),'XTickLabel',{[],'20%'})
        
        %% plotting the histogram across age groups per behavior
        if ag == length(agegroups)
            H3plot_beg = HMplot_beg+numSPcols;       H3plot_end = H3plot_beg+1;
            subplot(numSProws,numSPcols,H3plot_beg:H3plot_end)

            hist_ags = sum(sum(cat(3,array2disp{cls,:}),3),2)./(length(channelOrder)*sum(WTC_numTrls_agegroup(cls,:)));
            barh(log2(period_globalvar),hist_ags,'FaceColor','k')
            yvalues = 1./period_globalvar;
            Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
            set(gca,'YTick',log2(1./Yticks(:)))
            set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
                'YTickLabel',num2str(Yticks'),...
                'YDir','reverse', 'layer','top')
            set(gca,'XLim',[0 0.2],'Xtick',(0:0.2:0.2),'XTickLabel',{[],'20%'})
            set(gca,'FontSize',12,'FontWeight','bold')
        else
            continue
        end        
    
    end
    
    %% plotting the histogram across behaviors per age group
    if cls == size(array2disp,1)
        H2plot_beg = Hplot_idx+1;       H2plot_end = H2plot_beg+1;
        subplot(numSProws,numSPcols,H2plot_beg:H2plot_end)

        hist_clss = sum(sum(cat(3,array2disp{:,ag}),3),2)./(length(channelOrder)*sum(WTC_numTrls_agegroup(:,ag)));
        barh(log2(period_globalvar),hist_clss,'FaceColor','k')
        yvalues = 1./period_globalvar;
        Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
        set(gca,'YTick',[])
        set(gca,'YLim',log2([min(period_globalvar),max(period_globalvar)]),...
            'YDir','reverse', 'layer','top')
        set(gca,'XLim',[0 0.2],'Xtick',(0:0.2:0.2),'XTickLabel',{[],'20%'})
        set(gca,'FontSize',12,'FontWeight','bold')
    end
 
end

%% save the figure
enhanceFigure(myfig)
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'wtcohere_Signif_HighCohere_EEGleadingACC_acrosschans_V4_normed2SubjTrlSize_acrossbehaviors.tif';
% cd([serverPath1,'\Kinematics'])
print(myfig,'-dpng', filename, '-r600')

%% add color bar (separately)
fig_cb = figure;
colormap(flipud(colormap('gray')))
caxis([cmin cmax])
cbhdl = colorbar;
ylabel(cbhdl,'Occurences of Relevant Kinematic-to-Neural Activity Coherence')
% set(get(cbhdl,'ylabel'),'string','Occurences of Relevant Kinematic-to-Neural Activity Coherence','FontSize',12,'FontWeight','bold')
set(cbhdl,'Position',[0.023 0.354 0.015 0.279],'FontSize',12,'FontWeight','bold')
%%
filename_cb = 'wtcohere_Signif_HighCohere_EEGleadingACC_acrosschans_colorbar_V2.tif';
print(fig_cb,'-dpng', filename_cb, '-r600')
%}



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%
% Plotting only histograms across behaviors for each age group
myfig = figure(13);
agegroup_str = {'6 mos.','7-9 mos.','10-12 mos.','13-18 mos.','19-24 mos.'};
load('wtcperiodglobalvariable_v2.mat')
load('wtcbehaviorlistglobalvariable.mat')
frqbinsIDX = find(period_globalvar <= 1/4 & period_globalvar >= 1/10);
array2disp = WTC_cnts_agegroup;
for ag = 1:length(agegroups)
    %% plotting the histogram across behaviors per age group
    subplot(1,6,ag)
    numTrls_ag = sum(WTC_numTrls_agegroup(:,ag));
    numTrlChn_ag = length(channelOrder)*numTrls_ag;
    hist_clss = sum(sum(cat(3,array2disp{:,ag}),3),2)./numTrlChn_ag;
    bar(log2(period_globalvar),hist_clss,'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor',[0.5 0.5 0.5])
    hold('on')
    bar(log2(period_globalvar(frqbinsIDX)),hist_clss(frqbinsIDX),'FaceColor','k','BarWidth',1,'EdgeColor','k')
    grid('on')
    title([agegroup_str{ag},10,' [',num2str(numTrlChn_ag),' trials]'],...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold')
    xvalues = 1./period_globalvar;
    Xticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
    set(gca,'XTick',log2(1./Xticks(:)))
    xlim(log2([period_globalvar(1),period_globalvar(end-2)]))
    set(gca,'XDir','reverse', 'layer','top')
    set(gca,'XTickLabel',num2str(Xticks'))
    
    if ag==3
        xlabel('         Frequency (Hz)','FontName','Times New Roman',...
            'FontSize',14,'FontWeight','bold','HorizontalAlignment','left')
    end
    set(gca,'YLim',[0 0.1],'Ytick',(0:0.05:0.1),'YTickLabel',[])
    if ag == 1
        set(gca,'YLim',[0 0.1],'Ytick',(0:0.05:0.1),'YTickLabel',...
            (0:5:10))
        ylabel('% Trials w/Artifacts')
    end
    set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold')
    
%     % Add lines to emphasize frequencies of importance
%     histvalIDX = (0:0.05:0.1);
%     % for line at 4 Hz
%     xval1 = repmat(log2(1/4),length(histvalIDX),1);  
%     line(xval1,histvalIDX,'color','r','linestyle','--','linewidth',2)
%     % for line at 10 Hz
%     xval2 = repmat(log2(1/10),length(histvalIDX),1);  
%     line(xval2,histvalIDX,'color','r','linestyle','--','Linewidth',2) 
    
end

% add total histogram at the end
subplot(1,6,6)
numTrls_Total = sum(WTC_numTrls_agegroup(:));
numTrlChn_Total = length(channelOrder)*numTrls_Total;
hist_clss = sum(sum(cat(3,array2disp{:,:}),3),2)./numTrlChn_Total;
bar(log2(period_globalvar),hist_clss,'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor',[0.5 0.5 0.5])
hold('on')
bar(log2(period_globalvar(frqbinsIDX)),hist_clss(frqbinsIDX),'FaceColor','k','BarWidth',1,'EdgeColor','k')
grid('on')
title(['All Subjects',10,' [',num2str(numTrlChn_Total),' trials]'],...
    'FontName','Times New Roman','FontSize',14,'FontWeight','bold')
xvalues = 1./period_globalvar;
Xticks = fliplr(2.^(fix(log2(min(xvalues))):fix(log2(max(xvalues)))));
set(gca,'XTick',log2(1./Xticks(:)))
xlim(log2([period_globalvar(1),period_globalvar(end-2)]))
set(gca,'XDir','reverse', 'layer','top')
set(gca,'XTickLabel',num2str(Xticks'))

set(gca,'YLim',[0 0.1],'Ytick',(0:0.05:0.1),'YTickLabel',[])
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold')

% % Add lines to emphasize frequencies of importance
% histvalIDX = (0:0.05:0.1);
% % for line at 4 Hz
% xval1 = repmat(log2(1/4),length(histvalIDX),1);  
% line(xval1,histvalIDX,'color','r','linestyle','--','linewidth',2)
% % for line at 10 Hz
% xval2 = repmat(log2(1/10),length(histvalIDX),1);  
% line(xval2,histvalIDX,'color','r','linestyle','--','Linewidth',2)

set(myfig,'position',[100 550 1300 165]);
tightfig;

%% save the figure
cd('\\bmi-nas-01\Contreras-UH\Infantdata\Infantdata\code\Zachs_Infant_decoding_files\wtc-Zach-files\figures')
filename = 'wtcohere_Signif_HighCohere_EEGleadingACC_histogram_acrosstrialsANDchans_V4_phaseextendto180.tif';
% cd([serverPath1,'\Kinematics'])
print(myfig,'-dpng', filename, '-r300')
