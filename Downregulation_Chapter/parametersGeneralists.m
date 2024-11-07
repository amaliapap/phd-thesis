function p=parametersGeneralists(n)

p.delta  = 0.05;
p.alphaN = 0.972;
p.alphaL = 0.3;
p.alphaF = 0.018;
p.rNstar = 0.4;
p.rLstar = 7.5;
p.alphaJ = 1.5;
p.rho    = 0.4*1e-6;
p.cF = 30;
p.cLeakage = 0.03;
p.cR = 0.03;
p.rhoCN = 5.68;
p.gammaN=1;
p.gammaDOC=1;
p.epsilonL=0.8;
p.bL=0.08;
p.bF=0.3;
p.bDOC=0.3;
p.bN=0.3;
p.bg=0.2;
p.mMin=1.1623*1e-9;
p.mMax=1;
p.fTemp2=2;
p.fTemp15=1.5;
p.epsilonF=1;

deltax = (log(p.mMax) - log(p.mMin)) / n;
i = (1:n)';  % Create a column vector of indices
x = log(p.mMin) + (i - 0.5) * deltax;

p.m = exp(x);
p.r = (3 / (4 * pi) * p.m / p.rho).^(1/3);
p.nu = min(1,3 * p.delta ./ p.r);