% addpath('~/Documents/Source/NUMmodel/matlab')
%
n = [6 10 15 20];

incr=3;
% set colormaps of different shades for each size-group
cmapG=brewermap(length(n)+incr,'blues');
cmapD=brewermap(length(n)+incr,'greens');
cmapPC=brewermap(length(n)+incr,'OrRd');
% cmapPC=brewermap(length(n)+incr,'oranges');

% cmapAC=brewermap(length(n)+incr,'reds');
% cmapAC=brewermap(length(n)+incr,'purples');
cmapAC=brewermap(length(n)+incr,'RdPu');


% colors for productivity
cmap=flip(cmocean('deep',5));
%%
clf
NPP = zeros(1,4);
Btot = NPP;
ProdHTL = NPP;
height=12;
width=16;
set(figure,'Renderer','Painters','Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto')
tiledlayout(3,4,'TileSpacing','tight',Padding='compact');
set(gcf,"Color",'w')
set(groot,'defaultAxesFontSize',10)
nCopepods = 6;
nUni = 10;
nCAct = 3;
nCPass =2;
%
% Seasonal watercolumn number of unicellular groups:
%
nexttile([1,3]);

n = [6 10 15 20]; % number of unicellular groups

for i = 1:length(n)
    p = setupNUMmodel([0.2 5],logspace(0,3,nCAct),n(i),nCopepods,1);
    p = parametersWatercolumn(p);

    p.tEnd = 3*365;
    sim = calcFunctions( simulateWatercolumn(p,60,-15) );
    
    NPP(i) = mean(sim.ProdNet);
    Btot(i) = mean(sim.Bplankton);
    ProdHTL(i) = mean(sim.ProdHTL);
 
  
colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:), cmapPC(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:), cmapAC(i+1,:), cmapAC(i+1,:), p.colGroup{end}};

    ixTime = sim.t>365;
    for j = 1:p.nGroups
        B = squeeze( sum( sim.B .* reshape(sim.dznom,1,numel(sim.dznom),1),2) );
        ix = p.ixStart(j):p.ixEnd(j);
        loglog(sim.p.m(ix), mean( B(ixTime,ix-p.idxB+1),1) ./ log(p.mUpper(ix)./p.mLower(ix)),...
            'color',  colGroupFamilies{j},'linewidth',1.5)
        
        hold on
        drawnow
        ylim([1 1000])
        xlim([1e-8 2e4])
        %xlabel('Size ({\mu}g_C)')
        %ylabel('Mean Sheldon biomass ({\mu}g_C/m^2)')
    end
end
ylim([10 1000])
plotlabel('a',false);


nexttile
ax=gca;
set(gca,'YTickLabel',[])
yyaxis right
plot(n(2)*[1,1],[0.5,1.3],'k',LineStyle=':')
hold on
plot(n,Btot/Btot(2),'Color',cmapPC(3,:),'linewidth',2,LineStyle='-')
plot( n,NPP/NPP(2),'Color',cmap(4,:),'linewidth',2,LineStyle='-')
plot(n,ProdHTL/ProdHTL(2),'Color',cmap(2,:),'linewidth',2,LineStyle='-')

ylim([0.5 1.3])
ax=gca;
ax.YAxis(1).Color='k';
ax.YAxis(2).Color='k';

xlabel('No. uni. classes')

lgd=legend('','B_{tot}','NPP','Prod_{HTL}','box','off',Location='best');
lgd.ItemTokenSize(1)=10;
plotlabel('b',false);
axis tight
lgd.Location='southeast';
%%
% Seasonal watercolumn number of copepod stages:
%


nexttile([1,3]);

n = [4 6 10 15];
for i = 1:length(n)
    p = setupNUMmodel([0.2 5],logspace(0,3,nCAct),nUni,n(i),1);
    p = parametersWatercolumn(p);
    
    p.tEnd = 3*365;
    sim = calcFunctions( simulateWatercolumn(p,60,-15) );
    
    NPP(i) = mean(sim.ProdNet);
    Btot(i) = mean(sim.Bplankton);
    ProdHTL(i) = mean(sim.ProdHTL);
 colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:), cmapPC(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:), cmapAC(i+1,:), cmapAC(i+1,:), p.colGroup{end}};

    ixTime = sim.t>365;
    for j = 1:p.nGroups
        B = squeeze( sum( sim.B .* reshape(sim.dznom,1,numel(sim.dznom),1),2) );
        ix = p.ixStart(j):p.ixEnd(j);
        loglog(sim.p.m(ix), mean( B(ixTime,ix-p.idxB+1),1) ./ log(p.mUpper(ix)./p.mLower(ix)), 'color',  colGroupFamilies{j},'linewidth',1.5)
        
        hold on
        drawnow
        ylim([1 1000])
        xlim([1e-8 2e4])
        %xlabel('Size ({\mu}g_C)')
        ylabel('Mean Sheldon biomass ({\mu}g C m^{-2})')
    end
end
ylim([10 1000])
plotlabel('c',false);

nexttile
set(gca,'YTickLabel',[])
yyaxis right
plot(n(2)*[1,1],[0.5,1.3],'k',LineStyle=':')
hold on
plot(n,Btot/Btot(2),'Color',cmapPC(3,:),'linewidth',2,LineStyle='-')
plot( n,NPP/NPP(2),'Color',cmap(4,:),'linewidth',2,LineStyle='-')
plot(n,ProdHTL/ProdHTL(2),'Color',cmap(2,:),'linewidth',2,LineStyle='-')

ylim([0.5 1.3])
ax=gca;
ax.YAxis(1).Color='k';
ax.YAxis(2).Color='k';
xlabel('No. copepod stages')
hold on
plot(n(2)*[1,1],[0.5,1.3],'k:')
plotlabel('d',false);

%%
% Seasonal watercolumn number of copepod groups:
%
nexttile([1,3]);

n = [3 5 7 9];
for i = 1:length(n)
    switch n(i)
        case 3
            p = setupNUMmodel(...
            1,...
            logspace(1,3,nCAct-1),nUni,nCopepods,1);
            colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:), cmapAC(i+1,:), p.colGroup{end}};

        case 5
            p = setupNUMmodel(...
            logspace(log10(0.2),log10(5),nCPass),...
            logspace(0,3,nCAct),nUni,nCopepods,1);
            colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:),cmapPC(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:),cmapAC(i+1,:), cmapAC(i+1,:), p.colGroup{end}};

        case 7
            p = setupNUMmodel(...
            logspace(log10(0.2),log10(5),nCPass+1),...
            logspace(0,3,nCAct+1),nUni,nCopepods,1);
            colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:),cmapPC(i+1,:), cmapPC(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:),cmapAC(i+1,:),cmapAC(i+1,:), cmapAC(i+1,:), p.colGroup{end}};

        case 9
            p = setupNUMmodel(...
            logspace(log10(0.2),log10(5),nCPass+1),...
            logspace(0,3,2*nCAct),nUni,nCopepods,1);
            colGroupFamilies={cmapG(i+1,:), cmapD(i+1,:),cmapPC(i+1,:), cmapPC(i+1,:), cmapPC(i+1,:), cmapAC(i+1,:),cmapAC(i+1,:),cmapAC(i+1,:), cmapAC(i+1,:),cmapAC(i+1,:), cmapAC(i+1,:),  p.colGroup{end}};

    end
    p = parametersWatercolumn(p);

    % Scale HTL mortality:
    %setHTL(0.005 * (n/6)^0, 0.1, true, false);

    p.tEnd = 3*365;
    sim = calcFunctions( simulateWatercolumn(p,60,-15) );
    
    NPP(i) = mean(sim.ProdNet);
    Btot(i) = mean(sim.Bplankton);
    ProdHTL(i) = mean(sim.ProdHTL);
 
    ixTime = sim.t>365;
    for j = 1:p.nGroups
        B = squeeze( sum( sim.B .* reshape(sim.dznom,1,numel(sim.dznom),1),2) );
        ix = p.ixStart(j):p.ixEnd(j);
        loglog(sim.p.m(ix), mean( B(ixTime,ix-p.idxB+1),1) ./ log(p.mUpper(ix)./p.mLower(ix)), 'color',  colGroupFamilies{j},'linewidth',1.5)
        
        hold on
        drawnow
        ylim([1 1000])
        xlim([1e-8 2e4])
        xlabel('Size ({\mu}g C)')
        %ylabel('Mean Sheldon biomass ({\mu}g_C/m^2)')
    end
end
ylim([10 1000])
plotlabel('e',false);

nexttile
set(gca,'YTickLabel',[])
yyaxis right
plot(n(2)*[1,1],[0.5,1.3],'k',LineStyle=':')
hold on
plot(n,Btot/Btot(2),'Color',cmapPC(3,:),'linewidth',2,LineStyle='-')
plot( n,NPP/NPP(2),'Color',cmap(4,:),'linewidth',2,LineStyle='-')
plot(n,ProdHTL/ProdHTL(2),'Color',cmap(2,:),'linewidth',2,LineStyle='-')

ylim([0.5 1.3])
ax=gca;
ax.YAxis(1).Color='k';
ax.YAxis(2).Color='k';
xlabel('No. copepod groups')
plotlabel('f',false);
axis tight
%%
% setFigWidth(17)
% setFigHeight(13)
% 
% exportgraphics(gcf,'sensitivity.pdf')


    