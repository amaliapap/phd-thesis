function calcDerivativesGeneralists(this, u, dNdt, dDOCdt, dudt)
    this.mort2 = this.mort2constant * u;
    this.jPOM = (1 - remin2) * this.mort2 + (1 - reminF) * this.JCloss_feeding ./ this.m;

    for i = 1:this.n
        dNdt = dNdt + ((-this.JNreal(i) + (1 - this.f(i)) * this.JlossPassive(i) + this.JNlossLiebig(i) + reminF * this.JCloss_feeding(i)) / this.m(i) + remin2 * this.mort2(i)) * u(i) / rhoCN;

        dDOCdt = dDOCdt + ((-this.JDOCreal(i) + (1 - this.f(i)) * this.JlossPassive(i) + this.JCloss_photouptake(i) + reminF * this.JCloss_feeding(i)) / this.m(i) + remin2 * this.mort2(i)) * u(i);

        dudt(i) = (this.Jtot(i) / this.m(i) - this.mortpred(i) - this.mort2(i) - this.mortHTL(i)) * u(i);
    end
end
