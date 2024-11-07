% load simNov_all_in
idxTime=109:120;
idxDep=1:3;
jLavg=squeeze(mean(sim.jLrates(idxTime,:,:,idxDep,:),1));
Bavg=double(squeeze(mean(sim.B(idxTime,:,:,idxDep,:),1)));

jLtime=squeeze(sim.jLrates(idxTime,:,:,idxDep,:));
mortHTLtime=squeeze(sim.mortHTLrates(idxTime,:,:,idxDep,:));
Btime=squeeze(sim.B(idxTime,:,:,idxDep,:));
mass=log(sim.p.m(sim.p.idxB:end))';
sumjL=0;
sum_mortHTL=0;
for i=1:51 % Try weighing by the biomass also over depth
    sumjL= sumjL+squeeze(jLtime(:,:,:,idxDep,i).*Btime(:,:,:,idxDep,i)).*mass(i);
    sum_mortHTL= sum_mortHTL+squeeze(mortHTLtime(:,:,:,idxDep,i).*Btime(:,:,:,idxDep,i)).*mass(i);
end
sim.mNPP=squeeze(mean(exp(sumjL./squeeze(sum(jLtime.*Btime,5))),4)); % avergae over the top 170m
sim.mHTL=squeeze(mean(exp(sum_mortHTL./squeeze(sum(mortHTLtime.*Btime,5))),4));

% npp=squeeze(sim.ProdNetAnnual);
% nppTime=sim.ProdNet(idxTime,:,:);
% for iTime=1:12
%     mNPP_norm(iTime,:,:,:)=sum(nppTime(iTime,:,:).*mNPP(iTime,:,:,:),1)./sum(nppTime,1);
% end
% mNPP_norm_avg=squeeze(mean(mNPP_norm,1));
% mNPP_Avg=squeeze(mean(mNPP,1));