function [deltas,fracResp]=calcDeltaXgeneralist(jN,jL,jF,jDOC)
delta  = 0.05;
alphaN = 0.972;
alphaL = 0.3;
alphaF = 0.018;
rNstar = 0.4;
rLstar = 7.5;
alphaJ = 1.5;
rho    = 0.4*1e-6;
cF = 30;
cLeakage = 0.03;
cR = 0.03;
rhoCN = 5.68;
gammaN=1;
gammaDOC=1;
epsilonL=0.8;
bL=0.08;
bF=0.3;
bDOC=0.3;
bN=0.3;
bg=0.2;
mMin=1.1623*1e-9;
mMax=1;
fTemp2=2;
fTemp15=1.5;
% f = Q10^((T-Tref)/10.); % Tref=10, Q10=2
epsilonF=1;
this.epsilonF=epsilonF;
%%
this.n=1;

deltax = (log(mMax) - log(mMin)) / this.n;
i = (1:this.n)';  % Create a column vector of indices
x = log(mMin) + (i - 0.5) * deltax;

this.m = exp(x);

this.r = (3 / (4 * pi) * this.m / rho).^(1/3);
this.nu = 3 * delta ./ this.r;
this.Jresp = cR * alphaJ * this.m;
this.JlossPassive = cLeakage ./ this.r .* this.m;

for i = 1:this.n

    this.JN(i)   = jN(i).*this.m; % converto to ug N/day
    this.JDOC(i) = jDOC(i).*this.m; % converto to ug C/day
    this.JL(i)   = jL(i).*this.m;
    this.JF(i)   = jF(i).*this.m;

    this.Jmax = alphaJ * this.m .* (1.0 - this.nu);
    % this.flvl = this.epsilonF * this.AF*F /...% & ! Note: adding a small number in the
    % ((this.AF*F+eps) + fTemp2*this.JFmax) ;%  ! demonominator to avoid negative values if F = JFmax = 0.

    JmaxT = fTemp2 * this.Jmax(i);
    % this.JF = this.flvl * fTemp2*this.JFmax;
    Jnetp(i) = this.JL(i) * (1 - bL) + this.JDOC(i) * (1 - bDOC) + this.JF(i) * (1 - bF) - fTemp2 * this.Jresp(i);

    this.dN(i) = min(1, 1 / this.JN(i) * (Jnetp(i) - this.JF(i) * (1 + bg)) / (1 + bg + bN));

    if this.JL(i) > 0
        this.dL(i) = min(1, 1 / (this.JL(i) * (1 - bL)) * ((this.JN(i) + this.JF(i)) * (1 + bg) - this.JDOC(i) * (1 - bDOC) - this.JF(i) * (1 - bF) + fTemp2 * this.Jresp(i) + bN * this.JN(i)));
    else
        this.dL(i) = -1;
    end

    if this.dL(i) < 0
        this.dL(i) = 0;
        this.dDOC(i) = min(1, 1 / (this.JDOC(i) * (1 - bDOC)) * ((this.JN(i) + this.JF(i)) * (1 + bg) - this.JF(i) * (1 - bF) + bN * this.JN(i) + fTemp2 * this.Jresp(i)));
    else
        this.dDOC(i) = 1;
    end

    if this.dN(i) < 0
        this.dN(i) = 0;
        this.Jnet(i) = 1 / (1 + bg) * (this.dDOC(i) * this.JDOC(i) * (1 - bDOC) + this.dL(i) * this.JL(i) * (1 - bL) + this.JF(i) * (1 - bF) - fTemp2 * this.Jresp(i) - bN * this.dN(i) * this.JN(i));
        f = this.Jnet(i) / (this.Jnet(i) + JmaxT);
        this.JNlossLiebig(i) = (1 - f) * this.JF(i) - f * JmaxT;
    else
        this.Jnet(i) = 1 / (1 + bg) * (this.dDOC(i) * this.JDOC(i) * (1 - bDOC) + this.dL(i) * this.JL(i) * (1 - bL) + this.JF(i) * (1 - bF) - fTemp2 * this.Jresp(i) - bN * this.dN(i) * this.JN(i));
        f = this.Jnet(i) / (this.Jnet(i) + JmaxT);
        this.JNlossLiebig(i) = 0;
    end

    if (this.Jnet(i) + JmaxT) == 0
        f = 0;
    end

    this.JNreal(i) = this.dN(i) * (1 - f) * this.JN(i);
    this.JDOCreal(i) = this.dDOC(i) * (1 - f) * this.JDOC(i);
    this.JLreal(i) = this.dL(i) * (1 - f) * this.JL(i);
    this.JFreal(i) = (1 - f) * this.JF(i);
    this.Jtot(i) = f * JmaxT - (1 - f) * this.JlossPassive(i);

    this.JNtot(i) = this.JNreal(i) + this.JFreal(i);

    this.JCloss_feeding(i) = (1 - this.epsilonF) / this.epsilonF * this.JFreal(i);
    this.JCloss_photouptake(i) = (1 - epsilonL) / epsilonL * this.JLreal(i);
    this.Jresptot(i) = fTemp2 * this.Jresp(i) + bDOC * this.JDOCreal(i) + bL * this.JLreal(i) + bN * this.JNreal(i) + bF * this.JFreal(i) + bg * this.Jnet(i);

    this.f(i) = f;
  
end
fracResp=zeros(1,this.n);
if this.Jnet>0
    fracResp=this.Jresptot./this.Jnet;
end
this.jN = this.JNreal;
this.jDOC = this.JDOCreal;
this.JF = this.JFreal;
deltas=[mean(this.dN) mean(this.dL) mean(this.dDOC)];