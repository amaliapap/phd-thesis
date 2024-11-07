%
% Plot the annual average trophic level of higher level production
% in the upper layer of a global simulation
%
[~, idxG ]=find(sim.p.typeGroups==5); 
ixG  = (sim.p.ixStart( idxG (1)):sim.p.ixEnd( idxG (end))) -sim.p.idxB+1;
% diatoms
[~, idxD ]=find(sim.p.typeGroups==3); 
ixD  = (sim.p.ixStart( idxD (1)):sim.p.ixEnd( idxD (end))) -sim.p.idxB+1;

% passive copepods
[~, idxPC ]=find(sim.p.typeGroups==10); 
ixPC  = (sim.p.ixStart( idxPC (1)):sim.p.ixEnd( idxPC (end))) -sim.p.idxB+1;

% active copepods
[~,idxAC]=find(sim.p.typeGroups==11); % Find index of Active Copepods
ixAC = (sim.p.ixStart(idxAC(1)):sim.p.ixEnd(idxAC(end))) -sim.p.idxB+1;

ixGroups={ixG,ixD,ixPC,ixAC};


PHTL = zeros(length(sim.t),length(sim.x),length(sim.y));
theta=getTheta(sim.p);

nx = length(sim.x);
ny = length(sim.y);
nz = length(sim.z);

Flows_groups=zeros(length(sim.t),nx,ny,nz,length(ixGroups),length(ixGroups));
Eflow_pom=zeros(length(sim.t),nx,ny,nz,size(sim.B,5));
Eflow_htl=Eflow_pom;
Eflow_doc=Eflow_pom;

N = sim.N;
DOC = sim.DOC;
Si = sim.Si;
BB = sim.B;

for iTime = find(sim.t>sim.t(end)-365)
    for ix = 1:nx
        for iy = 1:ny
            for iz=1:nz
                if isfield(sim.p,'idxSi')
                    u = [squeeze(N(iTime,ix,iy,iz)), ...
                        squeeze(DOC(iTime,ix,iy,iz)), ...
                        squeeze(Si(iTime,ix,iy,iz)), ...
                        squeeze(BB(iTime,ix,iy,iz,:))'];
                else
                    u = [squeeze(N(iTime,ix,iy,iz)), ...
                        squeeze(DOC(iTime,ix,iy,iz)), ...
                        squeeze(BB(iTime,ix,iy,iz,:))'];
                end
                rates = getRates(sim.p, u, squeeze(sim.L(iTime,ix,iy,iz)), squeeze(sim.T(iTime,ix,iy,iz)));
                for i=1:sim.p.n-sim.p.idxB+1
                    Eflow_pom(iTime,ix,iy,iz,i)=rates.jPOM(i).*squeeze(BB(iTime,ix,iy,iz,i));
                    Eflow_htl(iTime,ix,iy,iz,i)= rates.mortHTL(i).*squeeze(BB(iTime,ix,iy,iz,i));
                    Eflow_doc(iTime,ix,iy,iz,i)= rates.jDOC(i).*squeeze(BB(iTime,ix,iy,iz,i));
                    for j=1:1:sim.p.n-sim.p.idxB+1
                        Eflow(i,j)= rates.jFreal(i)*theta(i,j)/sum(theta(i,:),2)*squeeze(BB(iTime,ix,iy,iz,i));
                    end
                end
                for iPred=1:length(ixGroups)
                    predGroup=cell2mat(ixGroups(iPred));
                    for jPrey=1:length(ixGroups)
                        preyGroup=cell2mat(ixGroups(jPrey));
                        Flows_groups(iTime,ix,iy,iz,iPred,jPrey)=sum(sum(Eflow(predGroup,preyGroup)));
                    end
                end
            end
        end
    end
end





