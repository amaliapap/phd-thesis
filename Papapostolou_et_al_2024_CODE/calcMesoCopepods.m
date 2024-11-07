function Bmeso_intTop170=calcMesoCopepods(sim)

p=sim.p;
ixC = (p.ixStart(p.typeGroups==10)-p.idxB+1):(p.ixEnd(p.typeGroups==11)-p.idxB+1);
rGroups=calcRadiusGroups(sim.p); % for copepods it returns length/2
r_mesoC=find(rGroups(ixC)>200/2 & rGroups(ixC)<=2e4/2);
Bmeso_C=squeeze(sum(sim.B(:,:,:,:,r_mesoC),5));
% Take depth-integrated in the top 170m
field = squeeze(sum(Bmeso_C(:,:,:,1:3),5));
dz = sim.dznom(1:3);
Bmeso_intTop170 =double(squeeze( sum(field.*reshape(dz ,1,1,1,numel(dz)),4) / 1000)); % g/m2
