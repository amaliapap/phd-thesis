%% Watercolumn and spectrum from global:

idx = calcGlobalWatercolumn(lat,lon,sim);
%%
iDepth=1;
simWC.ProdNetwc=sim.ProdNet(:,idx.x,idx.y);
%%
g_biom=51;
ProdHTLwcNewg=zeros(length(sim.t),g_biom);
simWC.p=sim.p;
simWC.t=sim.t;
time=1:length(simWC.t);
 for iTime = 1:length(time)
     i=iTime;
     simWC.B(i,:)   = squeeze(sim.B(iTime,idx.x, idx.y, iDepth,:));
     simWC.N(i)   = squeeze(sim.N(iTime,idx.x, idx.y, iDepth));
     simWC.Si(i)  = squeeze(sim.Si(iTime,idx.x, idx.y, iDepth));
     simWC.DOC(i) = squeeze(sim.DOC(iTime,idx.x, idx.y, iDepth));
     simWC.L(i) = squeeze(sim.L(iTime,idx.x, idx.y, iDepth));
     simWC.T(i) = squeeze(sim.T(iTime,idx.x, idx.y, iDepth));
     ixX = idx.x;
     ixY = idx.y;
     iDepth = 1;
     if isfield(sim,'Si')
         u = [sim.N(iTime,ixX, ixY,iDepth), ...
             sim.DOC(iTime,ixX, ixY,iDepth), ...
             sim.Si(iTime,ixX, ixY,iDepth), ...
             squeeze(sim.B(iTime,ixX, ixY, iDepth, :))'];
     else
         u = [sim.N(iTime,ixX, ixY,iDepth), ...
             sim.DOC(iTime,ixX, ixY,iDepth), ...
             squeeze(sim.B(iTime,ixX, ixY, iDepth, :))'];
     end
     rates = getRates(sim.p, u, sim.L(iTime,ixX,ixY,iDepth), sim.T(iTime,ixX,ixY,iDepth));
     simWC.rates.jLreal(i, :) = rates.jLreal;
     simWC.rates.jFreal(i, :) = rates.jFreal;
     simWC.rates.jDOC(i, :)   = rates.jDOC;
     simWC.rates.jN(i, :)     = rates.jN;
     simWC.rates.jSi(i, :)    = rates.jSi;
     simWC.rates.jTot(i, :)   = rates.jTot;
     simWC.rates.f(i, :)      = rates.f;
     simWC.rates.jF(i, :)     = rates.jF;
     simWC.rates.jMax(i, :)   = rates.jMax;
     simWC.rates.jRespTot(i, :) = rates.jRespTot;
     simWC.rates.jR(i, :)     = rates.jR;
     simWC.rates.jLossPassive(i, :) = rates.jLossPassive;
     simWC.rates.jPOM(i, :) = rates.jPOM;
     simWC.rates.mortpred(i, :) = rates.mortpred;
     simWC.rates.mortHTL(i, :) = rates.mortHTL;
     simWC.rates.mort2(i, :)   = rates.mort2;
     simWC.rates.mort(i, :)    = rates.mort;

 end
 %
 for iTime=1:length(time)
    ProdHTLwcNewg(iTime,:)=simWC.rates.mortHTL(iTime,:).*simWC.B(iTime,:);
    % NPPwcNew(iTime)=sum((simWC.rates.jLreal(i,:)-simWC.rates.jRespTot(i,:)).*simWC.B(i,:));
 end
 ProdHTLwcNew=sum(ProdHTLwcNewg,2);
 simWC.ProdHTLwc=ProdHTLwcNew;
 simWC.mHTLdepth =squeeze(mean(mHTL(time,idx.x,idx.y,:),1,'omitnan'));
 %% ChemostatTrophicEfficiencyPlots(simWC)
% [mass_rangeUni, CfixedUni,~,npp_og]   = CfixedRespTotalTime(simWC,simWC.rates,month);
simWC.ratesTime=extractRatesTime(simWC.rates,ixMonth);
% k=1;% depth
% iTime=month;
%     [~, ~,lambda_htl(i,k)]  = maxTrophicLevel(sim,rates,squeeze(sim.B( iTime,k,:)));
%             mNPP(i,k) = exp(sum(rates.jLreal.*squeeze(sim.B( iTime,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
%                 ./(sum(rates.jLreal.*squeeze(sim.B( iTime,k,:))))  );
%             mHTL(i,k) = exp(sum(rates.mortHTL.*squeeze(sim.B( iTime,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
%                 ./(sum(rates.mortHTL.*squeeze(sim.B( iTime,k,:))))  );
% WaterColumnTrophicEfficiencyPlots(simWC,month);
 cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')


