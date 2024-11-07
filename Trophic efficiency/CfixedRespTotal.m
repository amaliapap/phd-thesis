%
% Algorithm that interpolates carbon fixed by generalists and diatoms 
% and returns total carbon fixed at each size range of the spectrum (ugC/L/day)
% 
function [mass_rangeU, CfixedU,nppG,npp_og]=CfixedRespTotal(sim,rates,intPoints,lat,lon)
% watch out ln:20 and 27
arguments
    sim struct;
    rates struct;
    intPoints = 20;
    lat double = [];
    lon double = [];
end
p=sim.p;

% generalists
[~, idxG ] = find(sim.p.typeGroups==5);   % index of generalists
ixG        = (sim.p.ixStart( idxG (1)):sim.p.ixEnd( idxG (end)))-sim.p.idxB+1;
% mG         = p.m(ixG);                            % mass of generalists
Cfixed_G(ixG) = sim.rates.jLreal(ixG).*sim.B(end,ixG)';       % fixed carbon  
% Cfixed_G(ixG) = sim.rates.jLreal(ixG).*sim.B(end,ixG);       % fixed carbon  

% diatoms
[~, idxD]  = find(sim.p.typeGroups==3); 
ixD        = (sim.p.ixStart( idxD (1)):sim.p.ixEnd( idxD (end)))-sim.p.idxB+1;% index of diatoms
% mD         = p.m(ixD);   
Cfixed_G(ixD) = sim.rates.jLreal(ixD).*sim.B(end,ixD)';        % fixed carbon
% Cfixed_G(ixD) = sim.rates.jLreal(ixD).*sim.B(end,ixD);        % fixed carbon

% mUni       = [mG,mD];
jnet  = rates.f.*rates.jMax./(1-rates.f);
jnet(isnan(jnet))=0;
%    case 'Generalists'
    iGroup = 1;
        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');
        jResp(ixG) = rates.jR(ixG) + betaL*rates.jLreal(ixG) + ( betaN*rates.jN(ixG).*rates.jLreal(ixG)./(rates.jLreal(ixG)+rates.jDOC(ixG)))...
            + (betaG*jnet(ixG).*rates.jLreal(ixG)./(rates.jLreal(ixG)+rates.jDOC(ixG)+rates.jFreal(ixG)));
    % case 'Diatoms'
    iGroup = 2;
        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');
        jResp(ixD) = rates.jR(ixD) + betaL*rates.jLreal(ixD) + ( betaN*rates.jN(ixD).*rates.jLreal(ixD)./(rates.jLreal(ixD)+rates.jDOC(ixD)))...
            + (betaG*jnet(ixD).*rates.jLreal(ixD)./(rates.jLreal(ixD)+rates.jDOC(ixD)));
    
        jResp(isnan(jResp))=0;
    nppG = abs(Cfixed_G - jResp);
%% Step1: Interpolation 
% Augment mass and interpolate the corresponding C-fixed

nGroups     = 2;
nPoints     = nGroups*intPoints; %       # points for the interpolation
m_nPoints   = zeros(nGroups,nPoints); % initialize mass for these points/rows<-groups
CfixedGroup = m_nPoints;      % initialize Cfixed for the interrpolation
jResp_int   = m_nPoints;
npp_int     = m_nPoints;
for iGroup = 1:2            % check for only 2 groups
    ix  = p.ixStart(iGroup):p.ixEnd(iGroup);
    m   = p.m(ix);            % group mass
    ixB = ix-p.idxB+1;      % index of group's trophic level
    %  augmented mass-array for the different groups, discretized in nPoints
    m_nPoints(iGroup,:) = logspace(log10(min(m)), log10(max(m)), nPoints); %this may need to be cell-array 
                                                                           % when copepods are included
    val_to_interp  = Cfixed_G(ixB) ;
    val_to_interp2 = jResp(ixB).*sim.B(end,ixB) ;
    val_to_interp3 = nppG(ixB);
    % Interpolation
    val_interp  = (interp1(log(m), val_to_interp, log(m_nPoints(iGroup,:)), 'linear'));
    val_interp2 = (interp1(log(m), val_to_interp2, log(m_nPoints(iGroup,:)), 'linear'));
    val_interp3 = (interp1(log(m), val_to_interp3, log(m_nPoints(iGroup,:)), 'linear'));

    CfixedGroup(iGroup,:) = val_interp;
    jResp_int(iGroup,:)   = val_interp2;
    npp_int(iGroup,:)     = val_interp3;  
end

B = reshape(m_nPoints',[],1);   % reshape augmented mass-array to a vector
C = reshape(CfixedGroup',[],1); % reshape interpolated C-fixed to a vector
D = reshape(jResp_int',[],1);   % reshape interpolated jResp to a vector
E = reshape(npp_int',[],1);   % reshape interpolated nppU to a vector

%% Total carbon fixed in each mass range of unicellular plankton
%
mUni          = p.m(sim.p.ixStart(1):sim.p.ixEnd(2)); 
mass_rangeUni = logspace(log10(min(min(mUni))),log10(max(max(mUni))), nPoints);
CfixedU       = zeros(1,length(mass_rangeUni)-1);  % initialize fixed carbon at each size range
jRespU        = zeros(1,length(mass_rangeUni)-1); 
nppU          = zeros(1,length(mass_rangeUni)-1); 

bins_no       =zeros(1,length(mass_rangeUni)-1);
%C fixed for Unicellular plankton    
for ix =1:numel(m_nPoints)                 % iterate through all size classes for the 2 groups
    for iRange=1:length(mass_rangeUni)-1
        if B(ix)>=mass_rangeUni(iRange) && B(ix)<=mass_rangeUni(iRange+1)
            bins_no(iRange)= bins_no(iRange)+1;
            CfixedU(iRange) = (CfixedU(iRange)+C(ix));%/(mass_rangeUni(iRange+1)-mass_rangeUni(iRange)); % normalise for range width
            jRespU(iRange)  = (jRespU(iRange)+D(ix));%/(mass_rangeUni(iRange+1)-mass_rangeUni(iRange));
            nppU(iRange) = (nppU(iRange)+E(ix));
        end
    end
end
mass_rangeU=mass_rangeUni(1:end-1);
npp = (CfixedU-jRespU)./bins_no;
% npp_Uint = nppU./bins_no;
npp_og=sum(Cfixed_G-jResp.*sim.B(end,ixG(1):ixD(end)));

%% Diagnostic Figure
% figure(1)
% clf
% % tiledlayout(1,2)
%     m = p.m((1+p.idxB-1):end);
% 
% nexttile
% set(gcf,'color','w');
%     plot(m_nPoints(1,:),CfixedGroup(1,:),'r+')
% hold on
%     plot(m_nPoints(2,:),CfixedGroup(2,:),'b+')
%     plot(mass_rangeU,CfixedU,'go')
%     plot(mass_rangeU,CfixedU,'g')
%     plot(mass_rangeU,jRespU,'m')
%     plot(mass_rangeU,npp,'cd')
% xlabel('Mass (\mugC)')
% ylabel('Fixed carbon (\mugCl^{-1}day^{-1})')
% set(gca,'XScale','log')
% axis tight
% legend('C_{fixed.G}','C_{fixed.D}','C_{fixed.U}','','jResp\timesU','NPP')
% 
% % nexttile
% % set(gcf,'color','w');
% % hold on
% %     % plot(mass_rangeU,npp,'g--')
% %     % plot(mass_rangeU,npp_Uint,'r:')
% % 
% %     plot(m(1:10), npp_og(1:10),'m*')
% %     plot(m(ixD), npp_og(ixD),'b*')
% %     plot(m(ixG),npp_int_decimated(ixG),'ms')
% %     plot(m(ixD),npp_int_decimated(ixD),'bs')
% % 
% % xlabel('Mass (\mugC)')
% % ylabel('Fixed carbon (\mugCl^{-1}day^{-1})')
% % set(gca,'XScale','log')
% % axis tight
% % legend('NPP','NPP_{int}','NPP_{OG}','NPP_{decimated}')
% 
% % nexttile
% % set(gcf,'color','w');
% %     plot(mass_range2,total_mort,'--')
% % hold on
% %     plot(m,sim.rates.mortHTL.*sim.B(end,:)','o')
% % 
% % xlabel('Mass (\mugC)')
% % ylabel('Prodthtl  (\mugCl^{-1}day^{-1})')
% % set(gca,'XScale','log')
% % axis tight
% % legend('\Sigmamort-range','mort')
