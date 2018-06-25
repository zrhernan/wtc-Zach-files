function makeWTCfigure( Rsq, period, coi, wtcsig, Wxy, t, varargin )
%MAKEWTCFIGURE makes the same exact figure used within the wavelet 
%coherence function by Aslak Grinsted (see wtc.m).
% Separate plot function was generated in order to have multiple plots in
% one figure.
%   Adapted from Aslak Grinsted's 'wtc' script.
%   Revised by Zach Hernandez, University of Houston, 2016

% initialize sampling period
dt = unique(diff(t));
 n = length(t);

% in case time difference isn't exact
if length(dt) > 1; dt = dt(1); end

% in case t is not a column vector
if size(t,2) > 1; t=t'; end

%----------default arguments for the WTC Plot-----------
Args=struct('Dj',1/12, ...    % this will do 12 sub-octaves per octave
            'S0',2*dt,...    % this says start at a scale of twice 
            'J1',[],...
            'Mother','Morlet', ...
            'FreqAxis', 1,...  % default is to display in frequency, not period
            'MaxScale',[],...   %a more simple way to specify J1
            'ArrowDensity',[30 30],...
            'ArrowSize',1,...
            'ArrowHeadSize',1);
Args=parseArgs(varargin,Args,{'BlackandWhite'});

if isempty(Args.J1)
    if isempty(Args.MaxScale)
        Args.MaxScale=(n*.17)*2*dt; %auto maxscale
    end
    Args.J1=round(log2(Args.MaxScale/Args.S0)/Args.Dj);
end

ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
Args.ArrowHeadSize=Args.ArrowHeadSize*120/ad;


H=imagesc(t,log2(period),Rsq);%#ok
%[c,H]=safecontourf(t,log2(period),Rsq,[0:.05:1]);
%set(H,'linestyle','none')

set(gca,'clim',[0 1])

d2c = brewermap(100,'Blues');
colormap(d2c)
% HCB=colorbar;

% frequency or period axis?
if Args.FreqAxis == 1          % for frequency
    yvalues = 1./period;
    Yticks = fliplr(2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues)))));
    set(gca,'YTick',log2(1./Yticks(:)))
%     ylabel('Frequency (Hz)')
else                            % for period
    yvalues = period;
    Yticks = 2.^(fix(log2(min(yvalues))):fix(log2(max(yvalues))));
    set(gca,'YTick',log2(Yticks(:)))
    ylabel('Period') 
end

set(gca,'YLim',log2([min(period),max(period)]), ...
    'YDir','reverse', 'layer','top', ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top')

hold on

%phase plot
if ~isreal(Wxy)
    aWxy=angle(Wxy);
else
    aWxy = Wxy;
end
aaa=aWxy;
aaa(Rsq<.5)=NaN; %remove phase indication where Rsq is low
%[xx,yy]=meshgrid(t(1:5:end),log2(period));

phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);

if ~all(isnan(wtcsig))
    [c,h] = contour(t,log2(period),wtcsig,[1 1],'k');%#ok
    set(h,'linewidth',2)
end
%suptitle([sTitle ' coherence']);
tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
set(hcoi,'alphadatamapping','direct','facealpha',.5)
hold off
end % EOF

