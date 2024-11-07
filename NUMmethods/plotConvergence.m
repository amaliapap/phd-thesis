
x0=0;
idxS = calcGlobalWatercolumn(60,-15,sim);
idxE = calcGlobalWatercolumn(0,-15,sim);
idx=idxS;
N = squeeze(double(sim.N(:,idx.x, idx.y, idx.z)))';


y0=0;
height=12;
width=16;
fig=figure(18);
set(fig,'Renderer','Painters','Units','centimeters',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto','Name','Sensitivity stages & no. groups')
tiledlayout(2,2,'TileSpacing','tight',Padding='compact');
set(gcf,"Color",'w')
set(groot,'defaultAxesFontSize',10)

nexttile(1)
contourf(sim.t,-sim.z(idxS.z),log10(N),linspace(-2,3,20),'LineStyle','none')
