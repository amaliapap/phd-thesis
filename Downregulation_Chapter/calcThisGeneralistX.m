function this=calcThisGeneralistX(N,L,F,DOC,m,isPureAutotroph, withAffinities,specificRates)
arguments
    N 
    L 
    F 
    DOC 
    m = 1e-4;
    isPureAutotroph = false;
    withAffinities=false;
    specificRates=true;
end
% isPureAutotroph=true;
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
fTemp2=1;%2;
fTemp15=1;%1.5;
% f = Q10^((T-Tref)/10.); % Tref=10, Q10=2
epsilonF=1;
this.epsilonF=epsilonF;
%%
this.n=1;

% deltax = (log(mMax) - log(mMin)) / this.n;
% i = (1:this.n)';  % Create a column vector of indices
% x = log(mMin) + (i - 0.5) * deltax;

% this.m = exp(x);
this.m=m;
this.r = (3 / (4 * pi) * this.m / rho).^(1/3);
this.nu = min(1,3 * delta ./ this.r);
this.Jresp = cR * alphaJ * this.m;
this.JlossPassive = cLeakage ./ this.r .* this.m;
this.Jmax = alphaJ * this.m .* (1.0 - this.nu);

if withAffinities==true
    this.AN = alphaN * this.r.^(-2) ./ (1+(this.r./rNstar).^(-2)) .* this.m;
    this.AL = alphaL./this.r .* (1-exp(-this.r./rLstar)) .* this.m .* (1.d0-this.nu);
    this.AF = alphaF*this.m;
    this.JFmax = cF./this.r .* this.m;
end
for i = 1:this.n

    if withAffinities==false
        this.JN(i)   = N(i).*this.m;
        this.JDOC(i) = DOC(i).*this.m;
        this.JL(i)   = L(i).*this.m;
        this.JF(i)   = F(i).*this.m;
    else
        this.JN(i)   = gammaN * fTemp15 * this.AN(i)*N*rhoCN ;% Diffusive nutrient uptake in units of C/time
        this.JDOC(i) = gammaDOC * fTemp15 * this.AN(i)*DOC; % Diffusive DOC uptake, units of C/time
        this.JL(i)   = epsilonL * this.AL(i)*L ; % Photoharvesting

        this.flvl = this.epsilonF * this.AF*F /...% & % Note: adding a small number in the
            ((this.AF*F+eps) + fTemp2*this.JFmax) ;%  % demonominator to avoid negative values if F = JFmax = 0.

        this.JF = this.flvl * fTemp2*this.JFmax;
    end

    JmaxT = fTemp2 * this.Jmax(i);

    if isPureAutotroph==true
        this.JF=0;
    end
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
this.fracResp=zeros(1,this.n);
if this.Jnet>0
    this.fracResp=this.Jresptot./this.Jnet;
end
this.fracResp(this.fracResp>1)=1;
this.jN = this.JNreal;
this.jDOC = this.JDOCreal;
this.JF = this.JFreal;

if specificRates==true
    this.JN = this.JNreal./this.m;
    this.JDOC = this.JDOCreal./this.m;
    this.JFreal = this.JFreal./this.m;
    this.JLreal = this.JLreal./this.m;
    this.Jtot = this.Jtot./this.m;
    this.Jresptot = this.Jresptot./this.m;
    this.Jnet = this.Jnet./this.m;
    this.Jresp = this.Jresp./this.m;
    this.JlossPassive = this.JlossPassive./this.m;
    this.Jmax = this.Jmax./this.m;
end

% deltas=[mean(this.dN) mean(this.dL) mean(this.dDOC)];