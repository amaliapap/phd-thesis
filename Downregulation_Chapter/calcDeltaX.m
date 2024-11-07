%
% Input: specific flux rates jX (1/day)
%
function [deltas,fracResp]=calcDeltaX(jN,jL,jSi,jDOC)
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
epsilonL=0.8;
bL=0.08;
bSi=0.3;
bDOC=0.3;
bN=0.3;
bg=0.2;
mMin=1e-8;
mMax=1;
fTemp2=2;
fTemp15=1.5;
v=0.6;

%%
this.n=1;

deltax = (log(mMax) - log(mMin)) / this.n;
i = (1:this.n)';  % Create a column vector of indices
x = log(mMin) + (i - 0.5) * deltax;

this.m = exp(x);
this.r = (3 / (4 * pi) * this.m / rho/(1-v)).^(1/3);
this.nu = 6^(2/3) * pi^(1/3) * delta * (this.m / rho).^(-1/3) * ...
          (1 + v^(2/3)) / (1 - v)^(2/3);
%--------------------------------------- 
%    WATCH OUT!!! ONLY TESTING HERE!
% tmp=this.m(3);
% this.m=tmp;
%---------------------------------------
this.Jresp = cR * alphaJ * this.m;
this.JlossPassive = cLeakage ./ this.r .* this.m;

for i = 1:this.n

    this.JN(i)   = jN(i).*this.m;   % convert to ug N/day
    this.JDOC(i) = jDOC(i).*this.m; % convert to ug C/day
    this.JL(i)   = jL(i).*this.m;   % convert to ug C/day
    this.JSi(i)  = jSi(i).*this.m;  % convert to ug Si/day

    this.Jmax = alphaJ * this.m .* (1.0 - this.nu);

    JmaxT = fTemp2 * this.Jmax(i);
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
    Jnet(i) = max(0, 1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i))));
    Jlim(i) = Jnet(i);
end

if this.JL(i) > 0
    this.dL(i) = min(1, 1 / (this.JL(i) * (1 - bL)) * (Jlim(i) * (1 + bg + bSi + bN) - (1 - bDOC) * this.JDOC(i) + fTemp2 * this.Jresp(i)));
else
    this.dL(i) = -1;
end

if this.dL(i) < 0
    this.dL(i) = 0;
    this.dDOC(i) = min(1, 1 / (this.JDOC(i) * (1 - bDOC)) * (Jlim(i) * (1 + bg + bSi + bN) + fTemp2 * this.Jresp(i)));
end

Jnetp(i) = this.dL(i) * this.JL(i) * (1 - bL) + this.dDOC(i) * this.JDOC(i) * (1 - bDOC) - fTemp2 * this.Jresp(i);
Jnet(i) = max(0,1 / (1 + bg) * (Jnetp(i) - (bN * this.dN(i) * this.JN(i) + bSi * this.dSi(i) * this.JSi(i))));

if this.dSi(i) < 1 && this.dN(i) < 1
    Jlim(i) = Jnet(i);
    if this.JN(i) > 0
        this.dN(i) = Jlim(i) / this.JN(i);
    end
    if this.JSi(i) > 0
        this.dSi(i) = Jlim(i) / this.JSi(i);
    end
end
this.Jnet = Jnet;

f = (Jnet(i))/(Jnet(i)+ JmaxT);
if ((Jnet(i) + JmaxT)==0) 
   f=0;
end 
  

    this.JNreal(i) = this.dN(i) * (1 - f) * this.JN(i);
    this.JDOCreal(i) = this.dDOC(i) * (1 - f) * this.JDOC(i);
    this.JLreal(i) = this.dL(i) * (1 - f) * this.JL(i);
    this.JSireal(i) = (1 - f) * this.JSi(i);
    this.Jtot(i) =max(0, f * JmaxT - (1 - f) * this.JlossPassive(i));

    this.JCloss_photouptake(i) = (1 - epsilonL) / epsilonL * this.JLreal(i);
    this.Jresptot(i) = fTemp2 * this.Jresp(i) + bDOC * this.JDOCreal(i) + bL * this.JLreal(i) + bN * this.JNreal(i) + bSi * this.JSireal(i) + bg * this.Jnet(i);

    this.f(i) = f;
end
fracResp=zeros(1,this.n);
if Jnet>0
    fracResp=this.Jresptot./Jnet;
end
this.jN = this.JNreal;
this.jDOC = this.JDOCreal;
this.JSi = this.JSireal;
deltas=[mean(this.dN) mean(this.dL) mean(this.dSi) mean(this.dDOC)];