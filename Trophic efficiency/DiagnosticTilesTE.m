% function that returns 3 tiles  for a given water column
%
% input: sim,lat,lon,NPP_extracted
% functions called: plotWC, plotPanelSpectrum
function DiagnosticTilesTE(sim,idx,lat,lon,site,time,ProdNetwc,ProdHTLwcNew,ixMonth)

cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 
cmap2  = flip(cmocean('deep',100));
ccmap2=cmap2(2:end,:); 
iDepth=1; iTime=ixMonth; % needed for plotPanelSpectrum


t1=nexttile( );
plotWC(squeeze(sim.B(:,idx.x, idx.y, idx.z,:)),sim,sim.z(idx.z));

set(gca,'XTickLabel','');
colormap(ccmap2);
set(colorbar,'visible','off')
xlabel('')
if site>1
    ylabel('')
end
% ylabel(cbar, '\mug_C/l','FontSize',12)
maxDepth=max(sim.z(idx.z));
axis square
% ylim([-maxDepth 0]);%same depth as global wc
ylim([max(-2000,-maxDepth) 0]);

t1.TitleHorizontalAlignment = 'left';

cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

t2=nexttile( );
plot(time,ProdNetwc,'Color',cmap(4,:)) % extracted from WC
hold on
if site==1
    ylabel('Production (\mug_Cl^{-1}day^{-1})')
end
yyaxis  right
plot(time,ProdHTLwcNew,'Color',cmap(2,:)) % mannually calculated
xlabel('Time (days)')
axis tight
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(2,:);
MTEannual=mean(ProdHTLwcNew,1)/mean(ProdNetwc,1);
axis square
if site==1
    leg1=legend('NPP','Prod_{HTL}','box','off',Location='best');
    leg1.ItemTokenSize(1)=10;
end
t3=nexttile( );
box on
 % [~, iTime] = min(abs(sim.t-time));
set(gca,'XTickLabel','');
plotPanelSpectrum(sim,iTime,iDepth,lat,lon);
% ylim(ylimSpectrum)
xlabel('Mass (\mug C)')
if site>1
    ylabel('')
    if site==2
        tt3=title(t3,'h. ',FontWeight='normal');
    elseif site==3
        tt3=title(t3,'i. ',FontWeight='normal');
    end
end
legend('')
axis square

ax=gca;
ax.XTick=[10^-6 10^-2 100];
ax.XTickLabel={ '10^{-6}' '10^{-2}' '10^2'};
ax.YTick=[10^-4 10^-2 1];
ax.YTickLabel={ '10^{-4}' '10^{-2}' '10^0'};

str =num2str(round(MTEannual,3));
title_list={'a. Seasonally stratified','b. Upwelling','c. Oligotrophic'};
tt1=title(t1,string(title_list(site)),"FontWeight","normal");
if site==2
    tt1.Position(1)=80;
    title(t2,append('e. $\bar{\epsilon}_{\mu}$=',str),"FontWeight","normal",Interpreter="latex")
else
    tt1.Position(1)=50;
    title(t2,append('f. $\bar{\epsilon}_{\mu}$=',str),"FontWeight","normal",Interpreter="latex")
end
