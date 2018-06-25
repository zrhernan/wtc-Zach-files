function [ listofwtcfeats ] = detect_importantWTCsbyTrial( Rsq, period, coi, wtcsig, Wxy )
%% Generate binary list of frequencies match all criteria:
    % - are statistically significant (based on Monte Carlo tests)
    % - have squared magnitude coherence > 0.5
    % - is within the Cone of Influence
    % - phase is between 0 and 90 degrees (i.e. in-phase to EEG lagging acceleration), and
    % - covers 10% of the total time range (dependent on frequency due to the Cone of Influence
    
    listofwtcfeats = cell(size(Rsq,1),1);
    for cls = 1:size(Rsq,1)
        for chn = 1:size(Rsq,2)
            numTrls = length(Rsq{cls,chn});
            if numTrls==0,continue,end
            for trl = 1:numTrls
                % Region of Contours where coherence values are statistically significant
                signifRegions = double(wtcsig{cls,chn}{trl} >= 1);

                % Region of Cone of Influence
                outsidecoi = nan(size(wtcsig{cls,chn}{trl}));
                for s=1:length(period{cls,chn}{trl})
                    outsidecoi(s,:)=(period{cls,chn}{trl}(s)<=coi{cls,chn}{trl});
                end

                % Region of Phases where Phase is between 0 and -(+)180 
                aWxy=angle(Wxy{cls,chn}{trl});
                aaa=aWxy;
                aaa(Rsq{cls,chn}{trl} < 0.5) = NaN; %remove phase indication where Rsq is low
                leadingphaseRegions = aaa <= 0 & aaa > -pi | aaa == pi;

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
                
                if isempty(signif_10prcntoftime),
                    listofwtcfeats{cls}(:,chn,trl) = zeros(67,1);
                else
                    listofwtcfeats{cls}(:,chn,trl) = signif_10prcntoftime;
                end

            end
        end
    end
    
end     %EOF