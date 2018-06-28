function [ listofwtcfeats ] = detect_importantWTCsbyTrial_forCMC( Rsq, period, coi, wtcsig, Wxy, ts )
%% Generate binary list of frequencies match all criteria:
    % - are statistically significant (based on Monte Carlo tests)
    % - have squared magnitude coherence > 0.5
    % - is within the Cone of Influence
    % - phase is between 0 and -90 degrees (i.e. in-phase to EEG leading EMG), and
    % - covers between [-1,1] seconds around the movement onset
    
listofwtcfeats = nan(size(Rsq{1,1},1),size(Rsq,1),size(Rsq,2));
for chn = 1:size(Rsq,1)
    for trl = 1:size(Rsq,2)
        % Region of Contours where coherence values are 
        % statistically significant
        signifRegions = double(wtcsig{chn,trl} >= 1);

        % Region of Cone of Influence
        outsidecoi = nan(size(wtcsig{chn,trl}));
        for s=1:length(period{chn,trl})
            outsidecoi(s,:)=(period{chn,trl}(s)<=coi{chn,trl});
        end

        % Region of Phases where Phase is between 0 and 179 (phase 
        % direction reversed in display) and equal to -(+)180
        % degrees
        aWxy=angle(Wxy{chn,trl});
        aaa=aWxy;
        aaa(Rsq{chn,trl} < 0.5) = NaN; %remove phase indication where Rsq is low
        leadingphaseRegions = aaa >= 0 & aaa < pi | aaa == -pi;

        % Find the Intersection of all 3 regions
        combineRegions = signifRegions & outsidecoi & leadingphaseRegions;

        % Indicate periods (aka frequencies) where the region covers
       % [-0.5, 0.5] seconds pre/post movement onset 
        signif_1secprepostMO = zeros(size(combineRegions,1),1);
        time_thr = 5;
        ts1idx = find(ts == -1);        ts2idx = find(ts ==  1);
        for frq = 1:size(combineRegions,1)
            
            if sum(combineRegions(frq,ts1idx:ts2idx)) >= round(sum(outsidecoi(frq,ts1idx:ts2idx))*(time_thr/100))
                signif_1secprepostMO(frq) = 1;
            end
        end

        if isempty(signif_1secprepostMO)
            listofwtcfeats(:,chn,trl) = zeros(size(combineRegions,1),1);
        else
            listofwtcfeats(:,chn,trl) = signif_1secprepostMO;
        end

    end
end

    
end     %EOF