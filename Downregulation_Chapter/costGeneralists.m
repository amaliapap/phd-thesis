function this=costGeneralists(this)
bL=0.08;
bF=0.3;
bDOC=0.3;
bN=0.3;
bg=0.2;
fTemp2=1;%2;

% this.Jresptot(i) = fTemp2 * this.Jresp(i) + bDOC * this.JDOCreal(i) 
% + bL * this.JLreal(i) + bN * this.JNreal(i) + bF * this.JFreal(i) + bg * this.Jnet(i);

jLcost    = bL*this.JLreal;
jNcost    = bN*this.JN;
jFcost   = bF*this.JFreal;
jDOCcost  = bDOC*this.JDOC;
jNet_cost   = bg*this.Jnet; % I think this should be the right one
jXtot_costs = this.Jresptot;
%-------------------------------
% MB = Metabolic budget array
%-------------------------------
Ycosts =[jLcost;jNcost;jFcost;jDOCcost; jNet_cost; fTemp2*this.Jresp'];

% total_MC=this.Jresptot;
total_MC=sum(Ycosts,1);
this.MB = [this.Jnet; (this.JlossPassive'); total_MC];
% total_MB = sum(MB,1); % =

Ycosts =[jLcost;jNcost;jFcost;jDOCcost; jNet_cost; fTemp2*this.Jresp'];
this.MBdetails = [Ycosts; this.Jnet; this.JlossPassive'; ];