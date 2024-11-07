function seasonTiles(ratesSummer,p,lettersStr,legendVisible)
% figure(1)
% clf
% set(gcf,'color','w');
% set(groot,'defaultAxesFontSize',10)
% tiledlayout(8,1,'tilespacing','compact','padding','loose','TileIndexing', 'columnmajor')
%% ***************************************************************'
%---------------------------------------------------------
%                SUMMER TILES
%----------------------------------------------------------
cmap={'#3487c5','#FFBE00','#D90000','#44AB9A','#AB449A'};

ixU = (p.ixStart(1):p.ixEnd(2))-p.idxB+1;
ixC = (p.ixStart(3):p.ixEnd(p.nGroups-1))-p.idxB+1;

ymax_U=max(max(ratesSummer.jTot(ixU),max(ratesSummer.jDOC(ixU),ratesSummer.jLreal(ixU))) );
ymin_U=max(max(max(ratesSummer.jRespTot(ixU),ratesSummer.jLossPassive(ixU)),max(ratesSummer.mortpred(ixU),ratesSummer.mort2(ixU))));

ymax_C=ceil(max(ratesSummer.jF(ixC))*10)/10;
ymin_C=ceil(max(max(ratesSummer.jRespTot(ixC)),max(ratesSummer.mortHTL(ixC)))*10)/10;

m =p.m(p.idxB:end);
tline = linspace(m(1),1000, 100);

iGroup=1;
    t1=nexttile([2 1]);
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m =p.m(ix+p.idxB-1);
sLinetype=':'; % for losses
    %
    % Generalists
    %
    if ((p.typeGroups(iGroup)~=3) && (p.typeGroups(iGroup)~=4))
        cjTot=semilogx(m, ratesSummer.jTot(ix), 'color',[0.65 0.65 0.65], 'linewidth',3);
        hold on
        % semilogx(m, ratesSummer.jN(ix), 'b-','linewidth',1.5)
        cjF=semilogx(m, ratesSummer.jFreal(ix), 'color',cmap{3}, 'linewidth',1.5);
        cjDOC=semilogx(m, ratesSummer.jDOC(ix), 'color',cmap{5},'linewidth',1.5);
        semilogx(m, ratesSummer.jN(ix), 'b','linewidth',1.5)
        % cjmax=semilogx(m, ratesSummer.jMax(ix), 'k:');
        cjL=semilogx(m, ratesSummer.jLreal(ix), 'color',cmap{2},'linewidth',2);
        cjPred=semilogx(m, -ratesSummer.mortpred(ix), 'color',cmap{3},'linewidth',1.5, 'linestyle',sLinetype);
        cjRes=semilogx(m, -ratesSummer.jRespTot(ix), 'k', 'linewidth',1.5, 'linestyle',sLinetype);
        cjVL=semilogx(m, -ratesSummer.mort2(ix), 'b','linewidth',1.5, 'linestyle','-');
        cjPL=semilogx(m, -ratesSummer.jLossPassive(ix), 'color',[0 0.5 0],'linewidth',1.5, 'linestyle',sLinetype);
    %   loglog(m, ratesSummer.mortStarve(ix), 'b-o','linewidth',1.5)
        cjHTL=loglog(m, -ratesSummer.mortHTL(ix), 'color',cmap{4},'linewidth',1.5,'LineStyle','-');
        cLosstot=loglog(m, -ratesSummer.mortHTL(ix)-ratesSummer.mortpred(ix) -ratesSummer.jRespTot(ix)-ratesSummer.mort2(ix)-ratesSummer.jLossPassive(ix),...
            'color',[ 93, 109, 126, .9*250]/250,'linewidth',3);
        cjHTL.Marker="o";
        cjHTL.MarkerSize=1.7;
        cjVL.Marker="o";
        cjVL.MarkerSize=1.7;
        zeroLine=plot(tline,0*tline,'k:','LineWidth',.5);
      tt(iGroup)=title(append(lettersStr(iGroup),'. ',p.nameGroup{iGroup}),'FontWeight','normal',FontSize=10);

    end
        hold off
    axis tight
    xlim(calcXlim(p))
    ylim([-ymin_U*1.5 ymax_U*1.5])
    set(gca,'XTickLabel','')
    set(gca,'Box','off')
    %
    % Diatoms:
    %
    iGroup=2;
    ixd = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    md = p.m(ix+p.idxB-1);
         sLinetype = ":"; % For LOSSES
    %
    t2=nexttile([2 1]);
        semilogx(md, ratesSummer.jTot(ixd), 'color',[0.65 0.65 0.65], 'linewidth',3)
        hold on
        semilogx(md, ratesSummer.jDOC(ixd), 'color',cmap{5},'linewidth',1.5)
        cjN=semilogx(md, ratesSummer.jN(ixd),'b','linewidth',1.5);
        cjSi=semilogx(md, ratesSummer.jSi(ixd), 'color',cmap{4}, 'linewidth',1.5, 'linestyle','-.');
        semilogx(md, ratesSummer.jLreal(ixd), 'color',cmap{2},'linewidth',2)
        % semilogx(m, ratesSummer.jMax(ixd), 'k:')
         cjPred=semilogx(m, -ratesSummer.mortpred(ixd), 'color',cmap{3},'linewidth',1.5, 'linestyle',sLinetype);
         %semilogx(m, ratesSummer.jR(ixd), 'k-.', 'linewidth',1.5)
         cjRes=semilogx(m, -ratesSummer.jRespTot(ixd), 'k', 'linewidth',1.5, 'linestyle',sLinetype);
         cjVL=semilogx(m, -ratesSummer.mort2(ixd), 'b','linewidth',1.5, 'linestyle','-');
         cjPL=semilogx(m, -ratesSummer.jLossPassive(ixd), 'color',[0 0.5 0],'linewidth',1.5, 'linestyle',sLinetype);
         %    loglog(m, ratesSummer.mortStarve(ixd), 'b-o','linewidth',1.5)
         cjHTL=loglog(m, -ratesSummer.mortHTL(ixd), 'color',cmap{4},'linewidth',1.5,'LineStyle','-');
         cLosstot=loglog(m, -ratesSummer.mortHTL(ixd)-ratesSummer.mortpred(ixd) -ratesSummer.jRespTot(ixd)-ratesSummer.mort2(ixd)-ratesSummer.jLossPassive(ixd),...
            'color',[ 93, 109, 126, .9*250]/250,'linewidth',3);
         cjHTL.Marker="o";
         cjHTL.MarkerSize=1.7;
         cjVL.Marker="o";
         cjVL.MarkerSize=1.7;
         zeroLine=plot(tline,0*tline,'k:','LineWidth',.5);
    tt(iGroup)=title(append(lettersStr(iGroup),'. ',p.nameGroup{iGroup}),'FontWeight','normal',FontSize=10);

    hold off
    axis tight
    xlim(calcXlim(p))
    ylim([-ymin_U*1.5 ymax_U*1.5])
    set(gca,'XTickLabel','')
    set(gca,'Box','off')
%
%
%% -------------------------------- 
%      % Rates for COPEPODS
%--------------------------------
%
%
ttStr=strings(1,5);
i=1;
% for spring
    typGroup=p.typeGroups(3);
    t3=nexttile([2 1]);
for iGroup = 3:p.nGroups-1
    % if p.typeGroups(iGroup)>typGroup
        % tCsumm(i)=nexttile;
    % end
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m =p.m(ix+p.idxB-1);
    sLinetype=':';
    % ymax_C=max(ratesSummer.jF(ix) );
    % ymin_C=max(max(ratesSummer.jRespTot(ix)),max(ratesSummer.mortHTL(ix)));

    semilogx(m, ratesSummer.jF(ix), 'color',cmap{3}, 'linewidth',1.5)
    hold on
    % semilogx(m, p.jFmax(ix),'k:')
    semilogx(m, -ratesSummer.mortpred(ix), 'color',cmap{3},'linewidth',1.5, 'linestyle',sLinetype)
    semilogx(m, -ratesSummer.jRespTot(ix), 'k', 'linewidth',1.5, 'linestyle',sLinetype)
    cjHTL=semilogx(m, -ratesSummer.mortHTL(ix), 'color',cmap{4},'linewidth',1.5, 'linestyle','-');
    cjHTL.Marker="o";
    cjHTL.MarkerSize=1.7;
    zeroLine=plot(tline,0*tline,'k:','LineWidth',.5);
  tt(i)=title(append(lettersStr(iGroup),'. ',p.nameGroup{iGroup}(1:15)),'FontWeight','normal',FontSize=10);
    
    xlim(calcXlim(p))
    set(gca,'XTickLabel','')
    set(gca,'Box','off')
    ylim([-ymin_C ymax_C])
    yticks([0 ymax_C/2 ymax_C])
    yticklabels({0,num2str(ymax_C/2),num2str(ymax_C)})
    if p.typeGroups(iGroup+1)==11 && (p.typeGroups(iGroup)==10)
        t4=nexttile([2 1]);
        % tt=title(append(lettersStr(iGroup),'. ','Active copepods'),'FontWeight','normal');
         i=i+2;
    end
end
xticks([10^-8 10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 1 10 100 1000]);
xticklabels({'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2','10^3'});
% ylabel('Rates (day^{-1})','FontSize',10); % 3.31772e-10,1.873077066128071,-1
% 
%-------------------------------- LEGEND ---------------------------------
leg=legend([cjL,cjDOC,cjF,cjN,cjSi,cjTot,cjPL,cjPred,cjRes,cjVL,cjHTL,cLosstot,zeroLine],...
    {'Light','DOC','Feeding','N','Si','Division rate',...
    'Passive losses','Predation','Respiration','Viral lysis','HTL','Total losses','y=0'}, ...
    'box','off');
leg.FontSize=8;
leg.NumColumns=6;
leg.Location='layout';
leg.Position(1)=.1;%.1
leg.Position(2)=.05; %.05
leg.Visible=legendVisible;

