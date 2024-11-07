function plotWC(B,sim,z)
ylimit = [-500 0];
xlimit = [sim.t(end)-365 sim.t(end)];

z = [0; z];

B = squeeze(double(sum(B,3))); % summing all plankton groups
B(:,2:length(z)) = B; %because now length(z)=length(z)+1
% B(:,1) = B(:,1);

B(B < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
% contourf(sim.t,-z,log10(B'),[linspace(-2,3,10)],'LineStyle','none')

contourf(sim.t,-z,log10(B'),'LineStyle','none')
ylabel('Depth (m)')
%set(gca,'ColorScale','log')
%shading interp
axis tight
colorbar
%caxis([0.1 100])
colorbar('ticks',-2:3,'limits',[-2 3])

% ylim(ylimit)
xlim(xlimit)
clim([-2,3])
xlabel('Time (days)')
end