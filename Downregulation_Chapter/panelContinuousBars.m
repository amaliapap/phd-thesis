%
% Makes a panel of the various contributions to respiration for a group
% with beta-factors (generalists or diatoms).
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%
function panelContinuousBars(var_i,varL,combinedLeg)
% 

m=var_i;
% colormap
cmap=flip(cmocean('deep',7));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(4:end-1,:)];
% ccmap2=flip(ccmap21);

% or newcolors
palePink=[196, 120, 134]/250;
% cmap=brewermap(21,'BrBG');
% cmap3=cmap(13:18,:);
% cmap4=brewermap(10,'PRGn');
% ccmap=cmap4(3:6,:);
% ccmap21=[cmap3(1:end-2,:);[215, 219, 221]/250;ccmap];
% ccmap2=flip(ccmap21);


set(gca,'yscale','linear')
hold on


jPlot=varL.MBdetails./sum(varL.MBdetails,1); % 8x5
temp=jPlot(6,:);
jPlot(6,:)=jPlot(1,:); % replace basal with light cost (row position)
jPlot(1,:)=temp;
        basal=fillbetweenlines(m, 0*jPlot(1,:), jPlot(1,:), ccmap2(1,:));
        ncost=fillbetweenlines(m, jPlot(1,:), jPlot(1,:)+jPlot(2,:), ccmap2(4,:));
        siFeedcost=fillbetweenlines(m, sum(jPlot(1:2,:),1), sum(jPlot(1:3,:),1), ccmap2(6,:));
        doccost=fillbetweenlines(m, sum(jPlot(1:3,:),1), sum(jPlot(1:4,:),1), ccmap2(2,:));
        synthcost=fillbetweenlines(m, sum(jPlot(1:4,:),1), sum(jPlot(1:5,:),1), ccmap2(7,:));
        light=fillbetweenlines(m, sum(jPlot(1:5,:),1), sum(jPlot(1:6,:),1), ccmap2(3,:)); % synth cost
        netgro=fillbetweenlines(m, sum(jPlot(1:6,:),1), sum(jPlot(1:7,:),1), palePink);
        passlo=fillbetweenlines(m, sum(jPlot(1:7,:),1), sum(jPlot(1:8,:),1), ccmap2(8,:));

        ax=gca;
        % ax.XTick=15:60:365;
        % ax.XTickLabel={'J','M','M','J','S','N'};
        % xtickangle(360)
        if combinedLeg==true
            legend('jR','N cost','Si/F cost','DOC cost',...
    'Synth. cost','L cost','growth','Passive loss','Jtot');
            set(legend,'visible','off')
        end


xlim([min(m) max(m)])
ylim([0 1])

xlabel('Var')
ylabel('')

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
