%λ
%
% Calculate all diagnostic rates for a simulation
%
% If tDay<0, average over the last month of the simulation
% else calculates the average of the month (15 days before and after)
%
% In:
%    sim - simulation structure
%    tDay - time (-1 by default)
%    depth - depth layer (-1 by default)
%
% Out :
% %   Rates - jLreal, jFreal, jDOC, jN,jSi,jTot,f, jMax,jR 
%           jRespTot, jLossPassive, mortpred, mortHTL, mort2, mort 
%

function Rates = calcAllRatesAvg(sim,tDay,depth)

arguments
    sim struct
    tDay double = sim.t(end);
    depth double = 1;
end

    tStep=length(sim.t)/sim.p.tEnd;
    % tStep=1;
if tDay<0
    %Average over the last month
    if sim.t(end)<30 %if the simulation lasted less than one month then the average is done over all the simulation
        time = 1;
    else
        [~, time] = min(abs(sim.t-(sim.t(end)-29)));
    end
    time = time:length(sim.t);
elseif tDay==sim.t(end)
    [~, time] = min(abs(sim.t-tDay));
else % take the average of the month; 15 days before and after the given day
    % time=round((tDay-15)*tStep):round((tDay+15)*tStep);
    time= find(sim.t>tDay-15 & sim.t<tDay+15);
end

% if depth<0
%     depth=1:length(sim.z);
% end
%
%Initialisation
%
switch sim.p.nameModel
    case 'chemostat'
        jLreal = zeros(length(time),sim.p.n-sim.p.nNutrients);
    case 'watercolumn'
        jLreal = zeros(length(time),length(sim.z),sim.p.n-sim.p.nNutrients);
    case 'global'
        jLreal = zeros(length(time),length(sim.x),length(sim.y),length(sim.z),sim.p.n-sim.p.nNutrients);
end
jFreal = jLreal; jDOC = jLreal;
jN = jLreal;
jSi = jLreal;
jTot = jLreal;
f=jLreal;
jF = jLreal;
jMax = jLreal;
jR = jLreal;
jRespTot = jLreal;
jLossPassive = jLreal;
mortpred = jLreal;
mortHTL = jLreal;
mort2 = jLreal;
mort = jLreal;
%
% Calculation
%
for i = 1:length(time)

    iTime = time(i);

    switch sim.p.nameModel

        case 'chemostat'
            iDepth=depth;
            if isfield(sim,'Si')
                u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), sim.Si(iTime,iDepth), ...
                    squeeze(sim.B(iTime,:))];
            else
                u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), squeeze(sim.B(iTime,iDepth,:))'];
            end
            rates = getRates(sim.p, u, mean(sim.L), sim.T);
            jLreal(i,:) = rates.jLreal;
            jFreal(i,:) = rates.jFreal;
            jDOC(i,:) = rates.jDOC;
            jN(i,:) = rates.jN;
            jSi(i,:) = rates.jSi;
            jTot(i,:) = rates.jTot;
            f(i,:) = rates.f;
            jF(i,:) = rates.jF;
            jMax(i,:) = rates.jMax;
            jR(i,:) = rates.jR;
            jRespTot(i,:) = rates.jRespTot;
            jLossPassive(i,:) = rates.jLossPassive;
            mortpred(i,:) = rates.mortpred;
            mortHTL(i,:) = rates.mortHTL;
            mort2(i,:) = rates.mort2;
            mort(i,:) = rates.mort;
        case 'watercolumn'
            for iDepth = 1:length(sim.z)
                if isfield(sim,'Si')
                    u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), sim.Si(iTime,iDepth), ...
                        squeeze(sim.B(iTime,iDepth,:))'];
                else
                    u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), squeeze(sim.B(iTime,iDepth,:))'];
                end
                rates = getRates(sim.p, u, sim.L(iTime,iDepth), sim.T(iTime,iDepth));
                jLreal(i,iDepth,:) = rates.jLreal;
                jFreal(i,iDepth,:) = rates.jFreal;
                jDOC(i,iDepth,:) = rates.jDOC;
                jN(i,iDepth,:) = rates.jN;
                jSi(i,iDepth,:) = rates.jSi;
                jTot(i,iDepth,:) = rates.jTot;
                f(i,iDepth,:) = rates.f;
                jF(i,iDepth,:) = rates.jF;
                jMax(i,iDepth,:) = rates.jMax;
                jR(i,iDepth,:) = rates.jR;
                jRespTot(i,iDepth,:) = rates.jRespTot;
                jLossPassive(i,iDepth,:) = rates.jLossPassive;
                mortpred(i,iDepth,:) = rates.mortpred;
                mortHTL(i,iDepth,:) = rates.mortHTL;
                mort2(i,iDepth,:) = rates.mort2;
                mort(i,iDepth,:) = rates.mort;
            end
    
        case 'global'
            for ixX = 1:length(sim.x)
                for ixY = 1:length(sim.y)
                    for iDepth = 1:length(sim.z)
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
                        jLreal(i,ixX,ixY,iDepth,:) = rates.jLreal;
                        jFreal(i,ixX,ixY,iDepth,:) = rates.jFreal;
                        jDOC(i,ixX,ixY,iDepth,:) = rates.jDOC;
                        jN(i,ixX,ixY,iDepth,:) = rates.jN;
                        jSi(i,ixX,ixY,iDepth,:) = rates.jSi;
                        jTot(i,ixX,ixY,iDepth,:) = rates.jTot;
                        f(i,ixX,ixY,iDepth,:) = rates.f;
                        jF(i,ixX,ixY,iDepth,:) = rates.jF;
                        jMax(i,ixX,ixY,iDepth,:) = rates.jMax;
                        jR(i,ixX,ixY,iDepth,:) = rates.jR;
                        jRespTot(i,ixX,ixY,iDepth,:) = rates.jRespTot;
                        jLossPassive(i,ixX,ixY,iDepth,:) = rates.jLossPassive;
                        mortpred(i,ixX,ixY,iDepth,:) = rates.mortpred;
                        mortHTL(i,ixX,ixY,iDepth,:) = rates.mortHTL;
                        mort2(i,ixX,ixY,iDepth,:) = rates.mort2;
                        mort(i,ixX,ixY,iDepth,:) = rates.mort;
                    end
                end
            end
    end
end

% Averaging over time
Rates.jLreal = squeeze(mean(jLreal,1));
Rates.jFreal = squeeze(mean(jFreal,1));
Rates.jDOC = squeeze(mean(jDOC,1));
Rates.jN = squeeze(mean(jN,1));
Rates.jSi = squeeze(mean(jSi,1));
Rates.jTot = squeeze(mean(jTot,1));
Rates.f = squeeze(mean(f,1));
Rates.jF = squeeze(mean(jF,1));
Rates.jMax = squeeze(mean(jMax,1));
Rates.jR = squeeze(mean(jR,1));
Rates.jRespTot = squeeze(mean(jRespTot,1));
Rates.jLossPassive = squeeze(mean(jLossPassive,1));
Rates.mortpred = squeeze(mean(mortpred,1));
Rates.mortHTL = squeeze(mean(mortHTL,1));
Rates.mort2 = squeeze(mean(mort2,1));
Rates.mort = squeeze(mean(mort,1));

