cd 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab'

x0=1;
y0=2;
height=7;  
width=10;

fig_no=12;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Selectivity');
myFONT = 11;

clf
set(gcf,'color','w');
set(gca,'fontsize', myFONT)

myColor = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
mortHTL_new = [0.005];
mAdults = logspace(log10(0.2), log10(1000), 5);

% mixing_new = 0.0014;
mixing_new = 0.1;
mHTL_new =.1;
isvector=false;
if isvector==true
mHTL_new = logspace(-2,2,4);
end
for i = 1:length(mortHTL_new)
    % subplot(2,1,i)
    p = setupGeneric(mAdults);
    p = parametersChemostat(p);
    p.tEnd = 365*18;
    p.d = mixing_new;
    boolQua = false;
    boolDec = true;
    if isvector==true
        setHTL(mortHTL_new, mHTL_new(i), boolQua, boolDec);
    else
        setHTL(mortHTL_new, mHTL_new, boolQua, boolDec);
    end
    % sim = simulateChemostat(p,40);
    % p = sim.p;
    sLegend = {};
    ixTime = length(sim.t); % Defaults to last time step
    set(gca,'xscale','log')
    
    xlim(calcXlim(sim.p))
    hold off
    xlabel('Mass ({\mu}g C)', 'Fontsize', myFONT)
    % yyaxis right

    syms x

    XX = logspace(log10(1e-09), log10(max(sim.p.m)), 1000);
    if isvector==false

        fdec = (1 ./ (1+ (x./ mHTL_new).^(-2))) * (x./ mHTL_new).^(-0.25); % declining
        fquad = (1 ./ (1+ (x./ mHTL_new).^(-2)));%*(1./log(1/sim.p.mDelta)); % non-declining
        fdec_res = double(subs(fdec, x, XX));
        fquad_res = double(subs(fquad, x, XX));
        hold on
        semilogx(XX, fdec_res, '--k', 'LineWidth', 1)
        legendentries(sim.p.nGroups+2) = area(XX,fquad_res,'EdgeColor', 'none', 'FaceColor', myColor{i}, ...
            'ShowBaseLine', 'off', 'FaceAlpha', i*0.1);
    else
        fdec = (1 / (1+ (x/ mHTL_new(i))^(-2))) * (x/ mHTL_new(i))^(-0.25); % declining
        fquad = (1 / (1+ (x/ mHTL_new(i))^(-2)));%*(1./log(1/sim.p.mDelta)); % non-declining
        fdec_res(i,:) = double(subs(fdec, x, XX));
        fquad_res(i,:) = double(subs(fquad, x, XX));
        hold on
        semilogx(XX, fdec_res(i,:), '--k', 'LineWidth', 1)
        legendentries(sim.p.nGroups+2) = area(XX,fquad_res(i,:),'EdgeColor', 'none', 'FaceColor', myColor{i}, ...
            'ShowBaseLine', 'off', 'FaceAlpha', i*0.1);
    end
    
    ax = gca;
    ax.YColor = 'k';
    ylabel('HTL selectivity', 'Fontsize', myFONT)
    xlim([XX(1), XX(end)])
    ylim([0 1.5])
    sLegend{end+1} = 'HTL selectivity';
    % legend(legendentries, sLegend, 'location','northeast','box','off', 'NumColumns', 7)
end