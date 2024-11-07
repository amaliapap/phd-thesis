% Calculates downregulated uptake rates for diatom cells
% input: N,L,Si,DOC concentrations & 
% withAffinities (boolean) to select if the rates are calculated with
% affinities
% if withAffinities==false the input corresponds to jX fluxes
% Returns rates jX in units (1/day)
%
function this=calcThisDiatomsX(N,L,Si,DOC,m,withAffinities,specificRates)
arguments % ~withAffinities  (jN,jL,jSi,jDOC,false)
    N 
    L 
    Si 
    DOC 
    m = 1e-4;
    withAffinities= false; 
    specificRates=true;
end
delta  = 0.05;
alphaN = 0.972;
alphaL = 0.3;
rNstar = 0.4;
rLstar = 7.5;
alphaJ = 1.5;
rho    = 0.4*1e-6;
cLeakage = 0.03;
cR = 0.03;
rhoCN = 5.68;
rhoCSi = 3.4;
gammaN=1;
gammaDOC=1;
gammaSi=1;

epsilonL=0.8;
bL=0.08;
bSi=0.3;
bDOC=0.3;
bN=0.3;
bg=0.2;
mMin=1e-8;
mMax=1;
fTemp2=1;%2;
fTemp15=1;%1.5;
v=0.6;

%%
this.n=1;

% deltax = (log(mMax) - log(mMin)) / this.n;
% i = (1:this.n)';  % Create a column vector of indices
% x = log(mMin) + (i - 0.5) * deltax;
% 
% this.m = exp(x);
this.m=m;
this.r = (3 / (4 * pi) * this.m / rho/(1-v)).^(1/3);
this.nu =min(1, 6^(2/3) * pi^(1/3) * delta * (this.m / rho).^(-1/3) * ...
          (1 + v^(2/3)) / (1 - v)^(2/3));
if withAffinities==true
    this.AN = alphaN * this.r.^(-2) ./ (1+(this.r/rNstar).^(-2)).* this.m / (1-v);
    this.AL = alphaL./this.r.* (1-exp(-this.r/rLstar)) .* this.m .* (1-this.nu) / (1-v);
    this.AF = 0.d0;
    this.JFmax = 0.d0;
end
this.Jresp = cR * alphaJ * this.m;
this.JlossPassive = cLeakage ./ this.r .* this.m;
this.Jmax = alphaJ * this.m .* (1.0 - this.nu);

for i = 1:this.n

    if withAffinities==false
        this.JN(i)   = N(i).*this.m; % convert units to ugN/day
        this.JDOC(i) = DOC(i).*this.m;
        this.JL(i)   = L(i).*this.m;
        this.JSi(i)   = Si(i).*this.m;
    else
        this.JN(i) =  fTemp15* gammaN * this.AN(i)*N*rhoCN; %! Diffusive nutrient uptake in units of C/time
        this.JDOC(i) = fTemp15*gammaDOC * this.AN(i)*DOC; %! Diffusive DOC uptake, units of C/time
        this.JSi(i) = fTemp15* gammaSi * this.AN(i)*Si*rhoCSi;%! Diffusive Si uptake, units of C/time
        this.JL(i) =  epsilonL * this.AL(i)*L; % ! Photoharvesting
    end
    JmaxT = fTemp2*this.Jmax(i);  

    Jlim(i)=0.0;

    Jnetp(i) = this.JL(i) * (1 - bL) + this.JDOC(i) * (1 - bDOC)- fTemp2 * this.Jresp(i);

    % Calculation of down-regulation factors for N-uptake
if this.JN(i) > 0
    this.dN(i) = min(1, 1./this.JN(i) * Jnetp(i) / (1 + bg + bN + bSi));
else
    this.dN(i) = 1;
end

% Si-uptake
if this.JSi(i) > 0
    this.dSi(i) = min(1, 1./this.JSi(i) * Jnetp(i) / (1 + bg + bN + bSi));
else
    this.dSi(i) = 1;
end

this.dDOC(i) = 1;

% Further checks and calculations
if this.dN(i) == 1 && this.dSi(i) == 1
    if this.JN(i) < this.JSi(i)
        Jlim(i) = this.JN(i);
        if this.JSi(i) > 0
            this.dSi(i) = Jlim(i) / this.JSi(i);
        end
    elseif this.JSi(i) <= this.JN(i)
        Jlim(i) = this.JSi(i);
        if this.JN(i) > 0
            this.dN(i) = Jlim(i) / this.JN(i);
        end
    end
end

if this.dSi(i) == 1 && this.dN(i) < 1
    Jlim(i) = this.JSi(i);
    if this.JN(i) > 0
        this.dN(i) = Jlim(i) / this.JN(i);
    end
elseif this.dN(i) == 1 && this.dSi(i) < 1
    Jlim(i) = this.JN(i);
    if this.JSi(i) > 0
        this.dSi(i) = Jlim(i) / this.JSi(i);
    end
end

if this.dSi(i) < 1 && this.dN(i) < 1
    this.dL(i) = 1;
    Jnetp(i) = this.dL(i) * this.JL(i) * (1 - bL) + this.dDOC(i) * this.JDOC(i) * (1 - bDOC) - fTemp2 * this.Jresp(i);
    % Jnet(i) = max(0, 1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i))));
    Jnet(i) = 1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i)));
    Jlim(i) = Jnet(i);
end

if this.JL(i) > 0
    this.dL(i) = min(1, 1 / (this.JL(i) * (1 - bL)) * (Jlim(i) * (1 + bg + bSi + bN) - (1 - bDOC) * this.JDOC(i) + fTemp2 * this.Jresp(i)));
else
    this.dL(i) = -1;
end

if this.dL(i) < 0
    this.dL(i) = 0;
    this.dDOC(i) = max(0,min(1, 1 / (this.JDOC(i) * (1 - bDOC)) * (Jlim(i) * (1 + bg + bSi + bN) + fTemp2 * this.Jresp(i))));
end

Jnetp(i) = this.dL(i) * this.JL(i) * (1 - bL) + this.dDOC(i) * this.JDOC(i) * (1 - bDOC) - fTemp2 * this.Jresp(i);
% Jnet(i) = max(0,1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i))));
Jnet(i) = 1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i)));

if this.dSi(i) < 1 && this.dN(i) < 1
    Jlim(i) = Jnet(i);
    if this.JN(i) > 0
        this.dN(i) = max(0, Jlim(i) / this.JN(i));
    end
    if this.JSi(i) > 0
        this.dSi(i) = max(0,Jlim(i) / this.JSi(i));
    end
end
this.Jnetp = Jnetp;

this.Jnet = Jnet;

f = (Jnet(i))/(Jnet(i)+ JmaxT);
if ((Jnet(i) + JmaxT)==0) 
   f=0;
end 
  
    this.JNreal(i) = this.dN(i) * (1 - f) * this.JN(i);
    this.JDOCreal(i) = this.dDOC(i) * (1 - f) * this.JDOC(i);
    this.JLreal(i) = this.dL(i) * (1 - f) * this.JL(i);
    this.JSireal(i) = this.dSi*(i)*(1 - f) * this.JSi(i);
    % this.Jtot(i) =max(0, f * JmaxT - (1 - f) * this.JlossPassive(i));
    this.Jtot(i) = f * JmaxT - (1 - f) * this.JlossPassive(i);

    this.JCloss_photouptake(i) = (1 - epsilonL) / epsilonL * this.JLreal(i);
    this.Jresptot(i) = fTemp2 * this.Jresp(i) + bDOC * this.JDOCreal(i)...
    + bL * this.JLreal(i) + bN * this.JNreal(i) + bSi * this.JSireal(i) + bg * this.Jnet(i);

    this.f(i) = f;
end
this.fracResp=zeros(1,this.n);
if Jnet>0
    this.fracResp=this.Jresptot./Jnet;
end
this.fracResp(this.fracResp>1)=1;

this.JN = this.JNreal;
this.JDOC = this.JDOCreal;
this.JSi = this.JSireal;

if specificRates==true
    this.JN = this.JNreal./this.m;
    this.JDOC = this.JDOCreal./this.m;
    this.JSi = this.JSireal./this.m;
    this.JLreal = this.JLreal./this.m;
    this.Jtot = this.Jtot./this.m;
    this.Jresptot = this.Jresptot./this.m;
    this.Jnet = this.Jnet./this.m;
    this.Jresp = this.Jresp./this.m;
    this.JlossPassive = this.JlossPassive./this.m;
    this.Jmax = this.Jmax./this.m;
end
% deltas=[mean(this.dN) mean(this.dL) mean(this.dSi) mean(this.dDOC)];