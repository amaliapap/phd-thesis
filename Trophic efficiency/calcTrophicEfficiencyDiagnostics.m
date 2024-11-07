%
% Calculate the functions from any simulation. Only the last year is
% calculated for watercolumn and global simulations.
%
% In:
%   sim - simulation structure
%   bPrintSummary - prints out a summary
%
% Out:
%   sim - same as input but with fields of global function added:
%        mNPP : mean size of primary producers
%        mHTL : mean size ingested by higher trophic levels
%
function sim = calcTrophicEfficiencyDiagnostics(sim, options)

arguments
    sim struct;
    options.bPrintSummary = true;
end

%********* Watch out! hardcoded value!!!! **********
%  Initialize the length of mass_range
cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')

%********* What's out! hardcoded value!!!! **********

switch sim.p.nameModel

    case 'chemostat'
        if sim.p.nNutrients==3
            u = [sim.N(end), sim.DOC(end),sim.Si(end), sim.B(end,:)];
            [sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.ProdBact, sim.eHTL,...
                sim.Bpico, sim.Bnano, sim.Bmicro] = ...=
                getFunctions(u, sim.L, sim.T);
        else
            u = [sim.N(end), sim.DOC(end), sim.B(end,:)];
            [sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.ProdBact, sim.eHTL,...
                sim.Bpico, sim.Bnano, sim.Bmicro] = ...=
                getFunctions(u, sim.L, sim.T);
        end

        [~, ~,lambda_htl]     = maxTrophicLevel(sim,sim.rates(end,:),sim.B(end,:)');
        mNPP = exp(sum(sim.rates.jLreal'.*sim.B(end,:).*log(sim.p.m(sim.p.idxB:end)))...
            ./(sum(sim.rates.jLreal'.*sim.B(end,:)))  );
        mHTL = exp(sum(sim.rates.mortHTL'.*sim.B(end,:).*log(sim.p.m(sim.p.idxB:end)))...
            ./(sum(sim.rates.mortHTL'.*sim.B(end,:)))  );
        sim.mNPP=mNPP;
        sim.mHTL=mHTL;
        sim.lambda_htl=lambda_htl;
        
    case 'watercolumn'
        sLibName = loadNUMmodelLibrary();
        ixTime = find(sim.t>(max(sim.t)-365)); % Just do the last year
        nZ = length(sim.z);
        
        lambda     = zeros(length(ixTime),nZ,size(sim.B,3));
        lambda_htl = zeros(length(ixTime),nZ);
        mNPP       = lambda_htl; 
        mHTL       = lambda_htl;
        
        for iTime = ixTime
            i = iTime-ixTime(1)+1; % index of calculated values for last year
            % Integrate over depth:
            for k = 1:nZ
                if ~isnan(sim.N(iTime,k))
                    % Get values at each depth and time:
                    if sim.p.nNutrients==3
                        u = [squeeze(sim.N(iTime,k)), ...
                            squeeze(sim.DOC(iTime,k)), ...
                            squeeze(sim.Si(iTime,k)), ...
                            squeeze(sim.B(iTime,k,:))'];
                    else
                        u = [squeeze(sim.N(iTime,k)), ...
                            squeeze(sim.DOC(iTime,k)), ...
                            squeeze(sim.B(iTime,k,:))'];
                    end
                    rates = getRates(sim.p, u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    % [~, ~,lambda_htl(i,k)]  = maxTrophicLevel(sim,rates,squeeze(sim.B(iTime,k,:)));
                    [lambda(i,k,:),lambda_htl(i,k)]  = calcTrophicLevel(sim.p,rates,squeeze(sim.B(iTime,k,:)));

                    mNPP(i,k) = exp(sum(rates.jLreal.*squeeze(sim.B(iTime,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
                        ./(sum(rates.jLreal.*squeeze(sim.B(iTime,k,:))))  );
                    mHTL(i,k) = exp(sum(rates.mortHTL.*squeeze(sim.B(iTime,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
                        ./(sum(rates.mortHTL.*squeeze(sim.B(iTime,k,:))))  );
                end
            end
            % iTimenow = iTime - ixTime(1)+1;
        end
        sim.mNPP=mNPP;
        sim.mHTL=mHTL;
        sim.lambda=lambda;
        sim.lambda_htl=lambda_htl;
    case 'global'

        sLibName = loadNUMmodelLibrary();
        ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
        % Get grid volumes:
        load(sim.p.pathGrid,'dv','dz','dx','dy');
        ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells
       
        nX = length(sim.x);
        nY = length(sim.y);
        nZ = 1;%length(sim.z);
        [~, idxAC ]=find(sim.p.typeGroups==11);
        ixAC  = (sim.p.ixStart( idxAC (1)):sim.p.ixEnd( idxAC (end)));% index of Active Copepods
        lambda = zeros(length(1:12),nX,nY,nZ,sim.p.n-sim.p.idxB+1);
        lambda_htl = zeros(length(1:12),nX,nY,nZ);
        mNPP       = lambda_htl; 
        mHTL       = lambda_htl;
        % TLtime     = zeros(length(sim.t),nX,nY,nZ,sim.p.n-sim.p.idxB+1);
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

        for iTime = ixTime
            for i = 1:nX
                for j = 1:nY
                    for k = 1%:nZ only surface now
                        if isfield(p,'idxSi')
                            u = [squeeze(N(iTime,i,j,k)), ...
                                squeeze(DOC(iTime,i,j,k)), ...
                                squeeze(Si(iTime,i,j,k)), ...
                                squeeze(B(iTime,i,j,k,:))'];
                        else
                            u = [squeeze(N(iTime,i,j,k)), ...
                                squeeze(DOC(iTime,i,j,k)), ...
                                squeeze(B(iTime,i,j,k,:))'];
                        end
                        conv = squeeze(dz(i,j,k));
                        rates = getRates(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                        [lambda(iTime,i,j,k,:),lambda_htl(iTime,i,j,k)]  =...
                            calcTrophicLevel(sim.p,rates,squeeze(sim.B(iTime,i,j,k,:)));
                        
                        % [~,~,lambda_htl(iTime,i,j,k)] =...
                        %     maxTrophicLevel(sim,rates,squeeze(sim.B(iTime,i,j,k,:)));
                        % mNPP(iTime,i,j,k) =...
                        %     exp(sum(rates.jLreal.*squeeze(sim.B(iTime,i,j,k,:)).*log(sim.p.m(sim.p.idxB:end))')...
                        %     ./(sum(rates.jLreal.*squeeze( sim.B(iTime,i,j,k,:) ))));
                        % mHTL(iTime,i,j,k) =...
                        %     exp(sum(rates.mortHTL.*squeeze(sim.B(iTime,i,j,k,:))...
                        %     .*log(sim.p.m(sim.p.idxB:end))')...
                        %     ./(sum(rates.mortHTL.*squeeze(sim.B(iTime,i,j,k,:))))  );
                        %  Bvec  = squeeze(sim.B(iTime,i, j, k,:));
                        %  TLtime(iTime,i,j,k,:)=calcTrophicLevelamp(sim,rates,Bvec);
                    end
                end
            end
        end
end
sim.lambda_htl=lambda_htl;
sim.lambda=lambda;
% sim.mNPP=mNPP;
% sim.mHTL=mHTL;

% if mass_range_length==length(mass_range)
%     disp('mass length ok')
% else
%     disp('mass length wrong')
% end