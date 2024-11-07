%
% Calculates Global mats of ProdNet, ProdHTL, MTE throughout the year
% And plots time series of global totals
%

cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
sLibName = loadNUMmodelLibrary();
% Get grid volumes:
sim.p.pathGrid="C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\TMs\MITgcm_2.8deg\grid.mat";
load(sim.p.pathGrid,'dv','dz','dx','dy');
ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells
one_year_vec=1:length(sim.t);
% if length(sim.t)>12
%     one_year_vec= (length(sim.t)-12+1):length(sim.t);
% end
ProdNet = zeros(length(one_year_vec), length(sim.x), length(sim.y));
ProdHTL = ProdNet;

jLreal = zeros(length(one_year_vec),length(sim.x), length(sim.y), length(sim.z),sim.p.n-sim.p.idxB+1);
jF = zeros(length(one_year_vec),length(sim.x), length(sim.y), length(sim.z),sim.p.n-sim.p.idxB+1);
jDOC = zeros(length(one_year_vec),length(sim.x), length(sim.y), length(sim.z),sim.p.n-sim.p.idxB+1);

nTime = length(one_year_vec); % modify to only take final year
nX = length(sim.x);
nY = length(sim.y);
nZ = length(sim.z);

%
% Extract fields from sim:
%
N = sim.N;
DOC = sim.DOC;
if isfield(sim.p,'idxSi')
    Si = sim.Si;
else
    Si = 0;
end
B = sim.B;
L = sim.L;
T = sim.T;
p = sim.p;

for iTime =1: nTime
    for i = 1:nX
        for j = 1:nY
            ProdNettmp = 0;
            ProdHTLtmp = 0;
            for k = 1:nZ
                if isfield(p,'idxSi')
                    u = [squeeze(N(one_year_vec(iTime),i,j,k)), ...
                        squeeze(DOC(one_year_vec(iTime),i,j,k)), ...
                        squeeze(Si(one_year_vec(iTime),i,j,k)), ...
                        squeeze(B(one_year_vec(iTime),i,j,k,:))'];
                else
                    u = [squeeze(N(one_year_vec(iTime),i,j,k)), ...
                        squeeze(DOC(one_year_vec(iTime),i,j,k)), ...
                        squeeze(B(one_year_vec(iTime),i,j,k,:))'];
                end
                [~, ProdNet1,ProdHTL1,~, ~,~,~,~] = ...
                    getFunctions(u, L(one_year_vec(iTime),i,j,k), T(one_year_vec(iTime),i,j,k), sLibName);
                conv = squeeze(dz(i,j,k));
                ProdNettmp = ProdNettmp + ProdNet1*conv;
                ProdHTLtmp = ProdHTLtmp +ProdHTL1*conv;
                %eHTL = eHTL + eHTL1/length(sim.z);
            end
            ProdNet(iTime,i,j) = ProdNettmp;
            ProdHTL(iTime,i,j) = ProdHTLtmp;
        end
    end
end
sim.ProdNet = ProdNet;
sim.ProdHTL = ProdHTL;
sim.eHTL = sim.ProdHTL./sim.ProdNet;
%
% Global totals
%
calcTotal = @(u) sum(u(ix(:)).*dv(ix(:))); % mug/day
%%
for i = 1:length(one_year_vec)
    sim.ProdNetTotal(i) = sum(sum(squeeze(sim.ProdNet(i,:,:)).*dx.*dy)); % mgC/day
    sim.ProdHTLTotal(i) = sum(sum(squeeze(sim.ProdHTL(i,:,:)).*dx.*dy,'omitnan'));% mgC/day
end

%% Figures
cmap  = flip(cmocean('deep',5));
ccmap = cmap(2:end,:); 
%-----------------------
% figure specifications
%-----------------------
x0=0; %positions (no need to change)
y0=0;
width=16; %figure width in cm
height=6; %figure height in cm

fig=figure('Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
set(groot,'defaultLineLineWidth',2)
set(gcf,'color','w');
time=1:12;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor')
nexttile(1)
set(gca,"YScale","log")
    semilogy(1:nTime,sum(sum(ProdNet,2),3),'Color',cmap(4,:))
    hold on
    semilogy(1:nTime,sum(sum(ProdHTL,2,'omitnan'),3),'Color',cmap(2,:))
    ylabel('Production (\mugCl^{-1}day^{-1})')
    yyaxis  right
    plot(1:nTime,sum(sum(sim.ProdHTL,2),3)./sum(sum(sim.ProdNet,2),3,'omitnan'),':')
    xlabel('Time (months)')
    axis tight
    title('Global total',FontWeight='normal')
    axis square

%A(isnan(A))=0;
sim.ProdHTL(isnan(sim.ProdHTL))=0;
nexttile(2)
set(gca,"YScale","log")
    semilogy(1:nTime,mean(mean(sim.ProdNet,2),3),'Color',cmap(4,:))
    hold on
    semilogy(1:nTime,mean(mean(sim.ProdHTL,2),3),'Color',cmap(2,:))
    ylabel('Production (\mugCl^{-1}day^{-1})')
    yyaxis  right
    plot(1:nTime,mean(mean(sim.ProdHTL,2),3)./mean(mean(sim.ProdNet,2),3),':')
    xlabel('Time (months)')
    axis tight
    lgd=legend('NPP','Prod_{HTL}','\epsilon_{\mu}','box','off',Location='best',fontsize=10);
    lgd.ItemTokenSize(1)=10;
    title('Global average',FontWeight='normal')
    axis square