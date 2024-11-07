function this=costDiatoms(this)
bL=0.08;
bSi=0.3;
bDOC=0.3;
bN=0.3;
bg=0.2;
fTemp2=1;%2;


jLcost    = bL*this.JLreal;
jNcost    = bN*this.JN;
jSicost   = bSi*this.JSi;
jDOCcost  = bDOC*this.JDOC;
jNet_cost   = bg*this.Jnet; % I think this should be the right one
jXtot_costs = this.Jresptot;
%-------------------------------
% MB = Metabolic budget array
%-------------------------------
Ycosts =[jLcost;jNcost;jSicost;jDOCcost; jNet_cost; fTemp2*this.Jresp'];

% total_MC=this.Jresptot;
total_MC=sum(Ycosts,1);
% this.MB = [this.Jnet./this.m'; (this.JlossPassive')./this.m'; total_MC./this.m'];
this.MB = [this.Jnet; (this.JlossPassive'); total_MC];

% total_MB = sum(MB,1); % =

Ycosts =[jLcost;jNcost;jSicost;jDOCcost; jNet_cost; fTemp2*this.Jresp'];
% this.MBdetails = [Ycosts./this.m'; this.Jnet./this.m'; this.JlossPassive'./this.m'; ];
this.MBdetails = [Ycosts; this.Jnet; this.JlossPassive'; ];