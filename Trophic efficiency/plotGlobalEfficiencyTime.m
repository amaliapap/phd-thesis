%% Figures
cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 
%-----------------------
% figure specifications
%-----------------------
x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
height=6; %figure height in cm

fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
set(gcf,'color','w');
time=1:12;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor')
nexttile(1)
set(gca,"YScale","log")
    semilogy(1:nTime,sum(sum(ProdNet,2),3),'Color',cmap(4,:))
    hold on
    semilogy(1:nTime,sum(sum(ProdHTL,2,'omitnan'),3),'Color',cmap(2,:))
    ylabel('Production (\mugCl^{-1}day^{-1})')
    yyaxis  right
    plot(1:nTime,sum(sum(sim.ProdHTL,2),3)./sum(sum(sim.ProdNet,2),3,'omitnan'),':')
    xlabel('Time (months)')
    axis tight
    title('Global total',FontWeight='normal')
    axis square

%A(isnan(A))=0;
sim.ProdHTL(isnan(sim.ProdHTL))=0;
nexttile(2)
set(gca,"YScale","log")
    semilogy(1:nTime,mean(mean(sim.ProdNet,2),3),'Color',cmap(4,:))
    hold on
    semilogy(1:nTime,mean(mean(sim.ProdHTL,2),3),'Color',cmap(2,:))
    ylabel('Production (\mugCl^{-1}day^{-1})')
    yyaxis  right
    plot(1:nTime,mean(mean(sim.ProdHTL,2),3)./mean(mean(sim.ProdNet,2),3),':')
    xlabel('Time (months)')
    axis tight
    lgd=legend('NPP','Prod_{HTL}','\epsilon_{\mu}','box','off',Location='best',fontsize=10);
    lgd.ItemTokenSize(1)=10;
    title('Global average',FontWeight='normal')
    axis square