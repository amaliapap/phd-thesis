%
% Makes a panel of the various contributions to respiration for a group
% with beta-factors (generalists or diatoms).
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%
function panelRespirationA(p,rates, iGroup,combinedLeg)

ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m = p.m(ix+p.idxB-1);

% colormap
cmap=flip(cmocean('deep',7));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(4:end-1,:)];

jR = rates.jR(ix);
set(gca,'yscale','linear')
basal=fillbetweenlines(m, 0*jR, jR, ccmap2(1,:));
hold on

ymax = max(rates.jR);
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

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_F = betaF*rates.jFreal(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix); % because jR_g= betaG*jNet = betaG*f*jMax
%         jR_g_old = betaG*rates.jTot(ix); this is slightly smaller than
%         jR_g
        jTot = rates.jTot(ix);% rates.f(ix).*rates.jMax(ix);
        
        doc=fillbetweenlines(m, jR, jR+jR_DOC, ccmap2(2,:));
        light=fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, ccmap2(3,:));
        nutr=fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, ccmap2(4,:));
        feed=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_F, ccmap2(5,:));
        gro=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, ccmap2(7,:));
        totgro=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot,  [215, 219, 221]/250   );

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot);
        ax=gca;
        ax.XTick=[10^-8 10^-6 10^-4 10^-2 1];
        
        legend({'Basal','DOC','Light','Nutrients','Feeding','Synthesis costs','Division rate'});
        if combinedLeg
            set(legend,'visible','off')
        end

    case 'Diatoms'
        cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaSi = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bSi');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_Si = betaSi*rates.jSi(ix);
%         jR_g = betaG*rates.jTot(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix);
        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);

        doc=fillbetweenlines(m, jR, jR+jR_DOC, ccmap2(2,:));
        light=fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, ccmap2(3,:));
        nutr=fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, ccmap2(4,:));
        si=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_Si, ccmap2(6,:));
        gro=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, ccmap2(7,:));
        totgro=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+max(0,jTot),  [215, 219, 221]/250   );

        ymax = 1.05*max(jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+max(0,jTot));
        ax=gca;
        ax.XTick=[10^-8 10^-6 10^-4 10^-2 1];

        legend({'Basal','DOC','Light','Nutrients','Silicate','Synthesis costs','Division rate'});
        if combinedLeg
            set(legend,'visible','off')
            feed=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N, ccmap2(5,:));
        leg=legend([basal,doc,light,nutr,feed,si,gro,totgro],...
        {'Basal','DOC','Light','Nutrients','Feeding','Silicate','Synthesis costs','Division rate'}, ...
        'box','off');
        leg.Position=[0.681457746116643,0.696093496315487,0.302469135802469,0.254070556309362];
        end

end
cd(thisDir) % Ensure it returns to the directory with the plotting script

title(append('Respiration of ',lower(p.nameGroup{iGroup})),"FontWeight","normal")
%semilogx(m, rates.jTot,'k-','linewidth',2);
set(gca,'xscale','log')

xlim([min(p.m(p.idxB:end)) max(m)])
ylim([0 ymax])

xlabel('Cell mass ({\mu}g C)')
ylabel('Respiration (day^{-1})')

hold off
% %-------------------------------- LEGEND ---------------------------------
% if combinedLeg==true && iGroup==2
%     feed=fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N, ccmap2(4,:));
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
