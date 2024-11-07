 % case 'watercolumn'
function [mNPP,mHTL]=mNPPmHTLwc(sim,time)
arguments
    sim struct;
    time double = -1;
end

cd('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
sLibName = loadNUMmodelLibrary();
if time<0
    ixTime = find(sim.t>(max(sim.t)-365)); % Just do the last year
else
    ixTime=time;
end
% nZ = 1;%length(sim.z);

lambda_htl = zeros(length(ixTime),1);
mNPP       = lambda_htl;
mHTL       = lambda_htl;

for iTime = ixTime
    i = iTime ;%1:iTime-ixTime(1)+1; % index of calculated values for last year
    % Integrate over depth:
    % for k = 1:nZ
        if ~isnan(sim.N( iTime))
            % Get values at each depth and time:
            u = [squeeze(sim.N( iTime)), ...
                squeeze(sim.DOC( iTime)), ...
                squeeze(sim.Si( iTime)), ...
                squeeze(sim.B(iTime,:))];
            rates = getRates(sim.p, u, sim.L( iTime), sim.T( iTime), sLibName);
            [~, ~,lambda_htl(i)]  = maxTrophicLevel(sim,rates,squeeze(sim.B( iTime,:)));
            mNPP(i) = exp(sum(rates.jLreal.*squeeze(sim.B( iTime,:)').*log(sim.p.m(sim.p.idxB:end))')...
                ./(sum(rates.jLreal.*squeeze(sim.B( iTime,:)')))  );
            mHTL(i) = exp(sum(rates.mortHTL.*squeeze(sim.B( iTime,:)').*log(sim.p.m(sim.p.idxB:end))')...
                ./(sum(rates.mortHTL.*squeeze(sim.B( iTime,:)')))  );
        end
    % end
    % iTimenow = iTime - ixTime(1)+1;

end