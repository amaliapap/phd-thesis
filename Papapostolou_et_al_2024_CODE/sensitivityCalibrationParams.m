%

% Vary the three parameters that are used for the calibration

%

function sensitivityCalibrationParams(showFigsOnly)

if showFigsOnly==true
    load('data_for_sensitivity_figure.mat')
end
n = 8;

mortHTL = linspace(0.0,0.01,n);
kw = linspace(0.01,0.15,n);
u = linspace(5,25,n);

if showFigsOnly==false
    figure(1)
    [objHTL, errHTL] = varyHTL( mortHTL );
    [objkw, errkw] = varykw( kw );
    [obju, erru] = varyu( u );
end
%%

cmap=flip(cmocean('deep',5));
cmapPC=brewermap(7,'OrRd');
ccmap(1:3,:)=[cmap(4,:); cmap(2,:); cmapPC(3,:)];

x0=0;
y0=0;
width=16;
height=11;
figure_number = 2;           % Desired figure number
fig = figure(figure_number); % Create the figure with the specified number
set(fig, 'Renderer','Painters','Units','centimeters',...
    'Position',[x0 y0 width height],...
    'PaperPositionMode','auto','Name','Sensitivity calibration');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(gca,'defaultLineLineWidth',2)
figure(2)
clf
tiledlayout(2,3,TileSpacing="tight",Padding="tight",TileIndexing="columnmajor")
% HTL:
nexttile(1)
for i=1:3
    plot(mortHTL,objHTL(:,i),'Color',ccmap(i,:),'linewidth',2)
    hold on
end
plot(mortHTL,0*mortHTL,'k--')
ylim([-1.5 2])
leg1=legend({'Pico','POC','Copepods'});
leg1.ItemTokenSize(1)=10;
axis square
plotlabel('a',false);

nexttile(2)
for i=1:3
    plot(mortHTL,objHTL(:,i+3),'Color',ccmap(i,:),'linewidth',2)
    hold on
end
ylabel('NPP (mg C m^{-2}day^{-1})')
xlabel('\mu_{HTL,0} (L \mug C^{-3/4}day^{-1})')
leg2=legend({'Eutrophic','Oligotrophic','Seasonal'});
leg2.ItemTokenSize(1)=10;
plotlabel('d',false);
ylim([0 2000])
axis square

% kw:

nexttile(3)
for i=1:3
    plot(kw,objkw(:,i),'Color',ccmap(i,:), 'linewidth',2)
    hold on
end
plot(kw,0*kw,'k--')
xlim([min(kw) max(kw)])
ylim([-1.5 2])
axis square
plotlabel('b',false);

nexttile(4)
for i=1:3
    plot(kw,objkw(:,i+3),'Color',ccmap(i,:),'linewidth',2)
    hold on
end

xlabel('k_w (m^{-1})')
xlim([min(kw) max(kw)])
ylim([0 2000])
axis square
plotlabel('e',false);

%legend({'Eutrophic','Oligotrophic','Seasonal'})
% u
nexttile(5)
for i=1:3
    plot(u,obju(:,i),'Color',ccmap(i,:), 'linewidth',2)
    hold on
end
plot(u,0*u,'k--')
ylim([-1.5 2])
xlim([min(u) max(u)])
axis square
plotlabel('c',false);

nexttile(6)
for i=1:3
    plot(u,obju(:,i+3),'Color',ccmap(i,:),'linewidth',2)
    hold on
end
xlabel('Sinking velocity (m day^{-1})')
axis tight
ylim([0 2000])
axis square
plotlabel('f',false);

%%
    function [obj,err] = varyHTL(mortHTL)

        objExpected = [1 1 1 800 200 1000]; % Pico, POC, copepods; NPP , eutrophic oligotrophic and seasonal
        load NUMmodel.mat

        p = setupNUMmodel(bParallel= true); % A fast version of the NUM setup
        p = parametersGlobal(p,1);
        mHTL = .1;
        bHTLdecline = false;
        bHTLquadratic = true;
        simHTL = {};
        p.tEnd = 3*365;

        for i = 1:length(mortHTL)
            mortHTL(i)
            setHTL(mortHTL(i), mHTL, bHTLquadratic, bHTLdecline);
            simHTL{i} = simulateGlobal(p,sim,bCalcAnnualAverages=true, bVerbose=false);
            obj(i,:) = EvaluateRun(simHTL{i});
            err(i) = double(mean(abs(log(obj(i,:) ./ objExpected))));
            drawnow
            %siminit = simHTL{i};
        end
    end

    function [obj,err] = varykw(kw)

        objExpected = [1 1 1 800 200 1000]; % Pico, POC, copepods; NPP , eutrophic oligotrophic and seasonal
        load NUMmodel.mat

        p = setupNUMmodel(bParallel= true); % A fast version of the NUM setup
        p = parametersGlobal(p,1);
        sim_kw = {};
        p.tEnd = 3*365;

        for i = 1:length(kw)
            kw(i)
            p.kw = kw(i);
            sim_kw{i} = simulateGlobal(p,sim,bCalcAnnualAverages=true, bVerbose=false);
            obj(i,:) = EvaluateRun(sim_kw{i});
            err(i) = double(mean(abs(log(obj(i,:) ./ objExpected))));
            drawnow
            %siminit = simHTL{i};
        end
    end

    function [obj,err] = varyu(u)
        objExpected = [1 1 1 800 200 1000]; % Pico, POC, copepods; NPP , eutrophic oligotrophic and seasonal
        load NUMmodel.mat
        p = setupNUMmodel(bParallel= true); % A fast version of the NUM setup
        p = parametersGlobal(p,1);
        sim_u = {};
        p.tEnd = 3*365;

        for i = 1:length(u)
            u(i)
            setSinkingPOM(p,u(i));
            sim_u{i} = simulateGlobal(p,sim,bCalcAnnualAverages=true, bVerbose=true);
            obj(i,:) = EvaluateRun(sim_u{i});
            err(i) = double(mean(abs(log(obj(i,:) ./ objExpected))));
            drawnow
        end
    end

end