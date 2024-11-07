p = setupNUMmodel(bParallel=true);
p = parametersGlobal(p);
p.tEnd = 20*365;
% sim = simulateGlobal(p,bCalcAnnualAverages=true);
%%
x0=1;
y0=2;
height=10;  
width=16;
fig no=l;
fig=figure(fig no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Convergence');

clf
set(gcf,'color','w');
clf
tiledlayout(2,2,'TileSpacing','compact')
plot_N_and_B(60,-15,sim);
plot_N_and_B(0,-15,sim);

nexttile(1)
title('Nitrogen ({\mu}g N l^{-1})','fontweight','normal')
set(gca,'xticklabel','')
plotlabel('a',false);

nexttile(2)
title('Plankton biomass ({\mu}g C l^{-1})','fontweight','normal')
ylabel('')
set(gca,'yticklabel','')
set(gca,'xticklabel','')
xlabel('')
plotlabel('b',false);

nexttile(3)
xlabel('Time (years)')
plotlabel('c',false);

nexttile(4)
ylabel('')
set(gca,'yticklabel','')
xlabel('Time (years)')
plotlabel('d',false);

function plot_N_and_B(lat,lon,sim)
    idx = calcGlobalWatercolumn(lat,lon,sim);
    z = sim.z(idx.z);
    z = [0; z];
    
    nexttile
    N = squeeze(double(sim.N(:,idx.x, idx.y, idx.z)))';
    N = [N(1,:); N];
    N(N<=0) = 1e-8;
    contourf(sim.t/365,-z,N,logspace(-2,3,20),'LineStyle','none')
    set(gca, 'colorscale','log')
    clim(10.^[-2,3])
    ylabel('Depth (m)')
    ylim([-200 0])

    nexttile
    B = squeeze(sum(sim.B(:,idx.x,idx.y,idx.z,:),5));

    B(:,2:length(z)) = B;
    B(:,1) = B(:,2);
    B(B < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
    contourf(sim.t/365,-z,B',logspace(-2,3,20),'LineStyle','none')
    set(gca, 'colorscale','log')
    clim(10.^[-2,3])
    h = colorbar('ticks',10.^(-2:3));
    xlabel('Time (days)')
    ylim([-200 0])
end

setFigWidth(18)
setFigHeight(8)
exportgraphics(gcf,'convergence.pdf')