%
% Makes a panel of the various contributions to respiration for a group
% with beta-factors (generalists or diatoms).
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%
function panelRatesStacked(p,rates, iGroup,combinedLeg,ylab)
arguments
    p 
    rates 
    iGroup 
    combinedLeg 
    ylab = true
end
% ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
% m = p.m(ix+p.idxB-1);
ix = 1:365;

% colormap
cmap=flip(cmocean('deep',7));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(4:end-1,:)];

jR = rates.jResp(ix);
set(gca,'yscale','linear')
% basal=fillbetweenlines(ix, 0*jR, jR, ccmap2(1,:));
hold on

ymax = max(rates.jResp);
thisDir=pwd;
switch p.nameGroup{iGroup}
    case 'Generalists'
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaF = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bF');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');

        jR_L = rates.jLreal(ix);
        jR_N = rates.jN(ix);
        jR_DOC = rates.jDOC(ix);
        jR_F = rates.jFreal(ix);
        % jR_g = rates.f(ix).*rates.jMax(ix); % because jR_g= betaG*jNet = betaG*f*jMax
%         jR_g_old = betaG*rates.jTot(ix); this is slightly smaller than
%         jR_g
        jTot = rates.jTot(ix);% rates.f(ix).*rates.jMax(ix);
        
        doc=fillbetweenlines(ix, 0*jR_DOC, jR_DOC, ccmap2(2,:));
        light=fillbetweenlines(ix, jR_DOC, jR_DOC+jR_L, ccmap2(3,:));
        nutr=fillbetweenlines(ix, jR_DOC+jR_L, jR_DOC+jR_L+jR_N, ccmap2(4,:));
        feed=fillbetweenlines(ix, jR_DOC+jR_L+jR_N, jR_DOC+jR_L+jR_N+jR_F, ccmap2(5,:));
        % gro=fillbetweenlines(ix, jR_DOC+jR_L+jR_N+jR_F, jR_DOC+jR_L+jR_N+jR_F+jR_g, ccmap2(7,:));
        % totgro=fillbetweenlines(ix, jR_DOC+jR_L+jR_N+jR_F+jR_g, jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot,  [215, 219, 221]/250   );
        plot(ix,max(0,jTot),'Color',[215, 219, 221]/250,LineStyle='-.')

        ymax = max(jR_DOC+jR_L+jR_N);
        ax=gca;
        ax.XTick=15:60:365;
        ax.XTickLabel={'J','M','M','J','S','N'};
        xtickangle(360)
        if ylab==true
            ylabel('Rates (day^{-1})')
        end
        yyaxis right
        ax.YAxis(2).Color=[.5 .5 .5];
        
        legend({'DOC','Light','Nutrients','Feeding','Division rate'});
        if combinedLeg
            set(legend,'visible','off')
        end

    case 'Diatoms'
        % cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

        
        jR_L = rates.jLreal(ix);
        jR_N = rates.jN(ix);
        jR_DOC = rates.jDOC(ix);
        jR_Si = rates.jSi(ix);
%         jR_g = betaG*rates.jTot(ix);
        % jR_g = rates.f(ix).*rates.jMax(ix);
        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);

        doc=fillbetweenlines(ix, 0*jR_DOC, jR_DOC, ccmap2(2,:));
        light=fillbetweenlines(ix, jR_DOC, jR_DOC+jR_L, ccmap2(3,:));
        nutr=fillbetweenlines(ix, jR_DOC+jR_L, jR_DOC+jR_L+jR_N, ccmap2(4,:));
        si=fillbetweenlines(ix, jR_DOC+jR_L+jR_N, jR_DOC+jR_L+jR_N+jR_Si, ccmap2(6,:));
        % gro=fillbetweenlines(ix, jR_DOC+jR_L+jR_N+jR_Si, jR_DOC+jR_L+jR_N+jR_Si+jR_g, ccmap2(7,:));
        % totgro=fillbetweenlines(ix, jR_DOC+jR_L+jR_N+jR_Si+jR_g, jR_DOC+jR_L+jR_N+jR_Si+jR_g+max(0,jTot),  [215, 219, 221]/250   );
        plot(ix,max(0,jTot),'Color',[215, 219, 221]/250,LineStyle='-.')
        ymax = 1.05*max(jR_DOC+jR_L+jR_N+jR_Si);
        ax=gca;
        ax.XTick=15:60:365;
        ax.XTickLabel={'J','M','M','J','S','N'};
        xtickangle(360)  

        if ylab==true
            ylabel('Rates (day^{-1})')
        end
        yyaxis right
        ax.YAxis(2).Color=[.5 .5 .5];

        legend({'DOC','Light','Nitrogen','Silicate','Division rate'});
        if combinedLeg
            set(legend,'visible','off')
        end

end
cd(thisDir) % Ensure it returns to the directory with the plotting script

title(append(p.nameGroup{iGroup}),"FontWeight","normal")
%semilogx(ix, rates.jTot,'k-','linewidth',2);
% set(gca,'xscale','log')

xlim([min(ix) max(ix)])
ylim([0 ymax])

xlabel('Time (days)')
ylabel('Division rate (day^{-1})')

hold off
% %-------------------------------- LEGEND ---------------------------------
% if combinedLeg==true && iGroup==2
%     feed=fillbetweenlines(ix, jR_DOC+jR_L+jR_N, jR_DOC+jR_L+jR_N, ccmap2(4,:));
% 
%     leg=legend([basal,doc,nutr,light,feed,si,gro,totgro],...
%         {'Basal','DOC','Light','Nutrients','Feeding','Silicate','Synthesis costs','Division rate'}, ...
%         'box','off');
%     leg.FontSize=8;
%     leg.NumColumns=6;
%     leg.Location='layout';
% end

%------------------------------------------------------------------
%
% Fill between y1 (lower) and y2 (upper).
%
    function h = fillbetweenlines(x,y1,y2,color)

        % Grey is default color if no color argument given:
        if nargin==3
            color = 0.5*[1 1 1];
        end

        % possibly flip vectors:
        x = reshape(x,1,length(x));
        y1 = reshape(y1,1,length(y1));
        y2 = reshape(y2,1,length(y2));

        % If color is a scalar, make it into a grey-scale:
        if length(color)==1
            color = color*[1 1 1];
        end

        % Find lower limit:
        ymin = min( [ylim() y1] );
        if strcmp(get(gca, 'yscale'),'log') && (ymin<=0)
            ymin = min(y1(y1>0));
            y1(y1<=0) = ymin;
        end

        x = [x x(end:-1:1)];
        y = [y1 y2(end:-1:1)];

        h=fill(x,y,color);
        set(h,'edgecolor',color,'edgealpha',0);

    end
%     time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
% 
% saveas(gcf,['panelResp_',time,'.png'])
% p.rand
end
