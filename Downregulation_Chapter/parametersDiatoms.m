function p=parametersDiatoms(n)

p.delta  = 0.05;
p.alphaN = 0.972;
p.alphaL = 0.3;
p.rNstar = 0.4;
p.rLstar = 7.5;
p.alphaJ = 1.5;
p.rho    = 0.4*1e-6;
p.cLeakage = 0.03;
p.cR = 0.03;
p.rhoCN = 5.68;
p.rhoCSi = 3.4;
p.gammaN=1;
p.gammaDOC=1;
p.gammaSi=1;
p.epsilonL=0.8;
p.bL=0.08;
p.bSi=0.3;
p.bDOC=0.3;
p.bN=0.3;
p.bg=0.2;
p.mMin=1e-8;
p.mMax=1;
p.fTemp2=2;
p.fTemp15=1.5;
p.v=0.6;

deltax = (log(p.mMax) - log(p.mMin)) / n;
i = (1:n)';  % Create a column vector of indices
x = log(p.mMin) + (i - 0.5) * deltax;

p.m = exp(x);
p.r = (3 / (4 * pi) * p.m / p.rho/(1-p.v)).^(1/3);
p.nu =min(1, 6^(2/3) * pi^(1/3) * p.delta * (p.m / p.rho).^(-1/3) * ...
          (1 + p.v^(2/3)) / (1 - p.v)^(2/3));