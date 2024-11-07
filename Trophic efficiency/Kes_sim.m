p = setupNUMmodel(bParallel=true);
setHTL(0.05, 1 ,true, false); % "Quadratic" mortality; not declining
p = parametersGlobal(p);
p.tEnd=10*365;
%% Need rerun to get monthly output because the EvaluateRun script assumes
% monthly output:
% load('C:\Users\ampap\Downloads\Ken''s sim\NUMmodel paper figures\NUMmodel.mat')
% sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";
% for Linux
% sim.p.pathGrid="\home\amalia\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";


% sim = simulateGlobal(p,sim,bCalcAnnualAverages=true);
% sim = simulateGlobal(p,bCalcAnnualAverages=true);

sim.p=p;

%%

idxTime=109:120;
jLavg=squeeze(mean(sim.jLrates(idxTime,:,:,:,:),1));
Bavg=double(squeeze(mean(sim.B(idxTime,:,:,:,:),1)));

jLtime=squeeze(sim.jLrates(idxTime,:,:,:,:));
Btime=squeeze(sim.B(idxTime,:,:,:,:));
mass=log(sim.p.m(sim.p.idxB:end))';
sumjL=0;
for i=1:51
    sumjL= sumjL+squeeze(jLtime(:,:,:,:,i).*Btime(:,:,:,:,i)).*mass(i);
end
mNPP=exp(sumjL./sum(jLtime.*Btime,5));
npp=squeeze(sim.ProdNetAnnual);
nppTime=sim.ProdNet(idxTime,:,:);
for iTime=1:12
    mNPP_norm(iTime,:,:,:)=sum(nppTime(iTime,:,:).*mNPP(iTime,:,:,:),1)./sum(nppTime,1);
end
mNPP_norm_avg=squeeze(mean(mNPP_norm,1));
mNPP_Avg=squeeze(mean(mNPP,1));
%%
figure(1)
clf
set(gcf,'color','w');

tiledlayout(2,2)
nexttile
surface(sim.x,sim.y,squeeze(log10(mNPP_Avg(:,:,1))'))
colorbar
shading flat

nexttile
surface(sim.x,sim.y,npp')
colorbar
shading flat

nexttile
surface(sim.x,sim.y,squeeze(log10(mNPP_norm_avg(:,:,1))'))
colorbar
shading flat
%%
mNPP(iTime,i,j,k) =...
    exp(sum(rates.jLreal.*squeeze(sim.B(iTime,i,j,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
    ./(sum(rates.jLreal.*squeeze( sim.B(iTime,i,j,k,:) ))));
mHTL(iTime,i,j,k) =...
    exp(sum(rates.mortHTL.*squeeze(sim.B(iTime,i,j,k,:))...
    .*log(sim.p.m(sim.p.idxB:end))')...
    ./(sum(rates.mortHTL.*squeeze(sim.B(iTime,i,j,k,:))))  );
