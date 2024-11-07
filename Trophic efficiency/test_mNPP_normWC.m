jLtime=squeeze(sim.jLrates(idxTime,idx.x,idx.y,1,:));
Btime=squeeze(sim.B(idxTime,idx.x,idx.y,1,:));
mass=log(sim.p.m(sim.p.idxB:end))';
%%
sumjL=0;
for i=1:51
    sumjL= sumjL+squeeze(jLtime(:,i).*Btime(:,i)).*mass(i);
end
mNPP=exp(sumjL./sum(jLtime.*Btime,2));
npp=squeeze(sim.ProdNetAnnual);
nppTime=sim.ProdNet(idxTime,idx.x,idx.y);
for iTime=1:12
    mNPP_n(iTime)=sum(nppTime(iTime).*mNPP(iTime),1)./sum(nppTime,1);
end
mNPP_norm=mNPP_n';
mNPP_norm_avg=squeeze(mean(mNPP_norm,1));
mNPP_Avg=squeeze(mean(mNPP,1));


weig=nppTime/sum(nppTime);