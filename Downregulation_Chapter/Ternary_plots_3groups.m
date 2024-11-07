addpath("C:\Users\ampap\Downloads\github_repo(1)")

%% Case: Diatoms
% add_ternary_paths

wlimits = ternary_axes_limits(1,'l',0,'low','l',1,'high','r',0,'low', false);
%
% Rows of A,B,C triplets for uniformly-spaced data, assuming 31 grid
% points along each of the three axes, bounded by wlimits defined earlier
[jN,jL,jSi] = ternary_arrays( 4,wlimits );
%
jDOC=0.3*ones(1,10);
deltas=zeros(10,4);
dN=zeros(1,10);
dL=dN;
dSi=dN;
for i=1:10
    [deltas(i,:),frResp(i)]=calcDeltaX(jN(i),jL(i),jSi(i),jDOC(i)); % it is done for only 1 size-class
    dN(i)=deltas(i,1);
    dL(i)=deltas(i,2);
    dSi(i)=deltas(i,3);
end
frResp(frResp>1)=1;

%% -------------------------------
%       Figure for DIATOMS
%---------------------------------

x0=1; %positions (no need to change)
y0=2;
height=10; %figure height in cm
width=13;

fig_no=1;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','deltas diatoms');

clf
set(gcf,'color','w');
custom_labels={'j_N','j_L','j_{Si}'};

tiledlayout(2,2,"TileSpacing","loose")
t1=nexttile(1);
ternary_figure(jN,jL,dN',wlimits,custom_labels)
tt1=title('a. \delta_N',FontWeight='normal');
t1.TitleHorizontalAlignment = 'left'; 

t2=nexttile(2);
ternary_figure(jN,jL,dL',wlimits,custom_labels)
tt2=title('b. \delta_L', FontWeight='normal');
t2.TitleHorizontalAlignment = 'left'; 

t3=nexttile(3);
ternary_figure(jN,jL,dSi',wlimits,custom_labels)
tt3=title('c. \delta_{Si}','FontWeight','normal');
t3.TitleHorizontalAlignment = 'left'; 

t4=nexttile(4);
climit=[0 max(frResp)];
ternary_figure(jN,jL,frResp',wlimits,custom_labels,climit)
tt4=title('d. j_{Resp}/j_{Net}','FontWeight','normal');
t4.TitleHorizontalAlignment = 'left'; 

tt1.Position(2)=.95;
tt2.Position(2)=.95;
tt3.Position(2)=.95;
tt4.Position(2)=.95;

%% Case: Pure phototrophs
% add_ternary_paths
deltas=zeros(10,3);
dN=zeros(1,10);
dL=dN;
dDOC=dN;

wlimits = ternary_axes_limits(1,'l',0,'low','l',1,'high','r',0,'low', false);
%
% Rows of A,B,C triplets for uniformly-spaced data, assuming 31 grid
% points along each of the three axes, bounded by wlimits defined earlier
[jN,jL,jDOC] = ternary_arrays( 4,wlimits );
%
for i=1:10
    [deltas(i,:),frResp(i)]=calcDeltaXphototroph(jN(i),jL(i),jDOC(i)); % it is done for only 1 size-class
    dN(i)=deltas(i,1);
    dL(i)=deltas(i,2);
    dDOC(i)=deltas(i,3);
end
frResp(frResp>1)=1;
%% -------------------------------
%  Figure for PURE PHOTOTROPHS
%---------------------------------

x0=0; %positions (no need to change)
y0=0;
height=10; %figure height in cm
width=13;

fig_no=2;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','deltas phototroph');

custom_labels={'j_N','j_L','j_{DOC}'};

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontWeight','normal')
tiledlayout(2,2,TileSpacing="loose",Padding="loose")
t1=nexttile;
ternary_figure(jN,jL,dN',wlimits,custom_labels)
tt1=title('a. \delta_N','FontWeight','normal');
t1.TitleHorizontalAlignment = 'left'; 

t2=nexttile;
ternary_figure(jN,jL,dL',wlimits,custom_labels)
tt2=title('b. \delta_L','FontWeight','normal');
t2.TitleHorizontalAlignment = 'left'; 

t3=nexttile;
ternary_figure(jN,jL,dDOC',wlimits,custom_labels)
tt3=title('c. \delta_{DOC}',FontWeight='normal');
t3.TitleHorizontalAlignment = 'left'; 


t4=nexttile;
climit=[0 max(frResp)+.1];
ternary_figure(jN,jL,frResp',wlimits,custom_labels,climit)
tt4=title('d. j_{Resp}/j_{net}','FontWeight','normal');
% set(get(gca,'XLabel'),'Position',[.5 0 0])
t4.TitleHorizontalAlignment = 'left'; 

tt1.Position(2)=.95;
tt2.Position(2)=.95;
tt3.Position(2)=.95;
tt4.Position(2)=.95;

%% Case: Generalists
% add_ternary_paths

wlimits = ternary_axes_limits(1,'l',0,'low','l',1,'high','r',0,'low', false);
%
% Rows of A,B,C triplets for uniformly-spaced data, assuming 31 grid
% points along each of the three axes, bounded by wlimits defined earlier
[jN,jL,jF] = ternary_arrays( 4,wlimits );
%
jDOC=0.5*ones(1,10);
deltas=zeros(10,3);
dN=zeros(1,10);
dL=dN;
dDOC=dN;

for i=1:10
    [deltas(i,:),frResp(i)]=calcDeltaXgeneralist(jN(i),jL(i),jF(i),jDOC(i)); % it is done for only 1 size-class
    dN(i)=deltas(i,1);
    dL(i)=deltas(i,2);
    dDOC(i)=deltas(i,3);
end
frResp(frResp>1)=1;

%% -------------------------------
%     Figure for GENERALISTS
%---------------------------------

x0=1; %positions (no need to change)
y0=2;
height=10; %figure height in cm
width=13;

fig_no=3;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','deltas generalists');

clf
set(gcf,'color','w');

custom_labels={'j_N','j_L','j_F'};

tiledlayout(2,2)
t1=nexttile(1);
ternary_figure(jN,jL,dN',wlimits,custom_labels)
tt1=title('a. \delta_N',FontWeight='normal');
t1.TitleHorizontalAlignment = 'left'; 

t2=nexttile(2);
ternary_figure(jN,jL,dL',wlimits,custom_labels)
tt2=title('b. \delta_L', FontWeight='normal');
t2.TitleHorizontalAlignment = 'left'; 

t3=nexttile(3);
ternary_figure(jN,jL,dDOC',wlimits,custom_labels)
tt3=title('c. \delta_{DOC}','FontWeight','normal');
t3.TitleHorizontalAlignment = 'left'; 

t4=nexttile(4);
climit=[0 max(frResp)];
ternary_figure(jN,jL,frResp',wlimits,custom_labels,climit)
tt4=title('d. j_{Resp}/j_{Net}','FontWeight','normal');
t4.TitleHorizontalAlignment = 'left'; 

tt1.Position(2)=.95;
tt2.Position(2)=.95;
tt3.Position(2)=.95;
tt4.Position(2)=.95;