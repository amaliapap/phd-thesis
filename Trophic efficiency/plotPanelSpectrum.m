%
% Plot a size spectrum at a given time (day).
% If the simulation is a watercolumn then indicate also the depth layer.
% If the simulation is global then indicate depth layer, and latitude and longitude.
% 
function s = plotPanelSpectrum(sim, time, iDepth, lat, lon, showLeg)
arguments
    sim struct;
    time = sim.p.tEnd;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
    showLeg logical = true;
end
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
[~, iTime] = min(abs(sim.t-time));

switch sim.p.nameModel
    
    case 'chemostat'
        % Extract from a chemostat:
        s.B = sim.B;
        if isfield(sim,'Si')
            u = [sim.N(iTime), sim.DOC(iTime), sim.Si(iTime), squeeze(sim.B(iTime,:))];
        else
            u = [sim.N(iTime), sim.DOC(iTime), squeeze(sim.B(iTime,:))];
        end
        s.L = mean(sim.L);
        s.T = sim.T;
        
    case 'watercolumn'
        % Extract from a single water column:
        z = sim.z;
        s.B = squeeze(sim.B(:,iDepth,:));
        if isfield(sim,'Si')
            u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), sim.Si(iTime,iDepth), ...
                squeeze(s.B(iTime,:))];
        else
            u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), squeeze(s.B(iTime,:))];
        end
        s.L = sim.L(iTime,iDepth);
        s.T = sim.T(iTime,iDepth);
        
    case 'global'
        % Extract from global run:
        idx = calcGlobalWatercolumn(lat,lon,sim);
        z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
        s.B = squeeze(sim.B(:, idx.x, idx.y, iDepth, :));
        if isfield(sim,'Si')
            u = [sim.N(iTime,idx.x, idx.y,iDepth), ...
                sim.DOC(iTime,idx.x, idx.y,iDepth), ...
                sim.Si(iTime,idx.x, idx.y,iDepth), ...
                squeeze(sim.B(iTime,idx.x, idx.y, iDepth, :))'];
        else
            u = [sim.N(iTime,idx.x, idx.y,iDepth), ...
                sim.DOC(iTime,idx.x, idx.y,iDepth), ...
                squeeze(sim.B(iTime,idx.x, idx.y, iDepth, :))'];
        end
        s.L = sim.L(iTime,idx.x, idx.y, iDepth);
        s.T = sim.T(iTime,idx.x, idx.y, iDepth);
end

s.p = sim.p;
s.t = sim.t;
if ~isfield('sim','rates')
    sim.rates = getRates(sim.p, u, s.L, s.T);
end
s.rates = sim.rates;
%
% Setup tiles:
%
%
% Spectrum
%
panelSpectrum2(s,iTime,showLeg)
% xlabel('')
% set(gca,'XTickLabel','');

% if strcmp(sim.p.nameModel, 'watercolumn')
% 
% % sgtitle(['Day = ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176), ', depth: ', num2str(z(iDepth)), ' m']) 
% 
% end

%
% Plot Sheldon biomass spectrum. The biomasses are normalised by the log of
% the ratio between upper and lower masses in each bin
%
function panelSpectrum2(sim, ixTime,showLeg)
arguments
    sim struct;
    ixTime {mustBeInteger} = length(sim.t); % Defaults to last time step
    showLeg logical = true;
end

p = sim.p;
rates=sim.rates;
% %
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;
    %
    
    ixAve = find( sim.t > sim.t(end)/2 );
 
    set(gca,'xscale','log','yscale','log')
    hold on
end
%
% Community spectrum:
%
[mc, Bc] = calcCommunitySpectrum(sim.B, sim);
legendentries(1) = loglog(mc, Bc, 'linewidth', 4.5,'color',[0.7, 0.7, 0.7]);
sLegend{1} = 'Community spectrum';

%
% Group spectra:
%
for iGroup = 1:p.nGroups-1
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;
    %
    % Plot the spectrum:
    %
    sim.B(sim.B<=0) = 1e-100; % avoid negative values
    legendentries(iGroup+1) = ...
        loglog(m, exp( mean( log(sim.B(ixAve, ixB)./log(Delta)),1)), 'linewidth',2,...
        'color',p.colGroup{iGroup});
    sLegend{iGroup+1} = p.nameGroup{iGroup};
end
ylim([0.0001,500])
xlim(calcXlim(sim.p))
hold off
xticks=[10^-8 10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 1 10];

xlabel('Mass ({\mu}gC)')
ylabel('Biomass ({\mu}gCl^{-1})')% Sheldon biomass
%legendentries=[dum,legendentries];
%sLegend=[captionedstrat,sLegend];
if showLeg==true
lh = legend(legendentries, sLegend, 'Box','off');
% lh.Location='northeast';

lh.ItemTokenSize(1) = 8;
lh.Location='southoutside';
% lh.Location='best';

% lh.Position = [0.071797526626239,0.035340684878997,0.291838837010124,0.221487597099021];
lh.FontSize = 8;
% 0.428049769932749,0.094644508313384,0.140104163282861,0.206832865708863
lh.NumColumns=4;
end