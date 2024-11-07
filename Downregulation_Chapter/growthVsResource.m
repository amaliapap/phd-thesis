addpath('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab')
n=10; % No. of size classes
incr=3;
% set colormaps of different shades for each size-group
cmapG=brewermap(n+incr,'blues');
cmapD=brewermap(n+incr,'greens');
cmapPC=brewermap(n+incr,'OrRd');

cmapAC=brewermap(n+incr,'RdPu');
nN=150;
nL=100;
nDOC=50;
nF=10;
nSi=150;
%% ----------------------------------
% Change resources for GENERALISTS
%------------------------------------
% Changing N
N=linspace(0,50,10);
for i=1:length(N)
    this{i} = generalistsf90(N(i),nL,nDOC,nF);
    jTot_mean(i)=mean(this{i}.Jtot./this{i}.m');
    jTot_all(i,:)=this{i}.Jtot./this{i}.m';
end

% Changing  DOC
DOC=linspace(0,300,10);
for i=1:length(DOC)
    this{i} = generalistsf90(nN,nL,DOC(i),nF);
    jTot_mean2(i)=mean(this{i}.Jtot./this{i}.m');
    jTot_all2(i,:)=this{i}.Jtot./this{i}.m';
end

% Changing  L
L=linspace(0,300,10);
for i=1:length(L)
    this{i} = generalistsf90(nN,L(i),nDOC,nF);
    jTot_mean3(i)=mean(this{i}.Jtot./this{i}.m');
    jTot_all3(i,:)=this{i}.Jtot./this{i}.m';
end
% Changing  F
F=linspace(0,100,10);
for i=1:length(L)
    this{i} = generalistsf90(nN,nL,nDOC,F(i));
    jTot_mean4(i)=mean(this{i}.Jtot./this{i}.m');
    jTot_all4(i,:)=this{i}.Jtot./this{i}.m';
end

%% ----------------------------------
%  Change resources for DIATOMS
%------------------------------------
% Changing N
N=linspace(0,50,10);
for i=1:length(N)
    thisD{i} = diatomsf90(N(i),nL,nDOC,nSi);
    jTot_meanD(i)=mean(thisD{i}.Jtot./thisD{i}.m');
    jTot_allD(i,:)=thisD{i}.Jtot./this{i}.m';
end

% Changing  DOC
DOC=linspace(0,300,10);
for i=1:length(DOC)
    thisD{i} = diatomsf90(nN,nL,DOC(i),nSi);
    jTot_meanD2(i)=mean(thisD{i}.Jtot./thisD{i}.m');
    jTot_allD2(i,:)=thisD{i}.Jtot./thisD{i}.m';
end

% Changing  L
L=linspace(0,300,10);
for i=1:length(L)
    thisD{i} = diatomsf90(nN,L(i),nDOC,nSi);
    jTot_meanD3(i)=mean(thisD{i}.Jtot./thisD{i}.m');
    jTot_allD3(i,:)=thisD{i}.Jtot./thisD{i}.m';
end
% Changing  Si
Si=linspace(10,150,10);
for i=1:length(L)
    thisD{i} = diatomsf90(nN,nL,nDOC,Si(i));
    jTot_meanD4(i)=mean(thisD{i}.Jtot./thisD{i}.m');
    jTot_allD4(i,:)=thisD{i}.Jtot./thisD{i}.m';
end
%%
cmap={'#3487c5','#FFBE00','#D90000','#44AB9A','#AB449A'};
cmap=flip(cmocean('deep',5));
cmapPC=brewermap(7,'OrRd');
ccmap(1:3,:)=[cmap(4,:); cmap(2,:); cmapPC(3,:)];

x0=1; %positions (no need to change)
y0=2;
height=11; %figure height in cm
width=14;

%***********************************
%   Generalists figure
%***********************************
fig_no=1;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','growth rates');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

t1=nexttile();
j_mean=plot(N, jTot_mean,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(N, jTot_all(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('Ν (\mug N l^{-1})')
ylabel('jTot')
ylim([0 2.5])

leg1=legend(j_mean,'$\overline{j_{Tot}}$',interpreter='latex',location='best',fontsize=12);
leg1.ItemTokenSize(1)=10;
colormap(cmapAC)
axis tight
axis square

t2=nexttile();
j_mean=plot(DOC, jTot_mean2,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
    plot(DOC, jTot_all2(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('DOC (\mug C l^{-1})')
ylabel('jTot')
ylim([0 2.5])

cbar=colorbar;
colormap(cmapAC)
ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square

t3=nexttile();
j_mean=plot(L, jTot_mean3,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(L, jTot_all3(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('L (\muE)')
ylabel('jTot')
ylim([0 2.5])

colormap(cmapAC)
ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square

t4=nexttile();
j_mean=plot(F, jTot_mean4,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(F, jTot_all4(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('F (\mug C l^{-1})')
ylabel('jTot')
cbar=colorbar;
colormap(cmapAC)

ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square
ylim([0 2.5])

sgtitle('Generalists')
%% 
%***********************************
%        Diatoms figure
%***********************************
fig_no=2;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','growth rates');

clf
set(gcf,'color','w');
set(groot,'defaultAxesFontSize',10)
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

t1=nexttile();
j_meanD=plot(N, jTot_meanD,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(N, jTot_allD(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('Ν (\mug N l^{-1})')
ylabel('jTot')

leg1=legend(j_meanD,'$\overline{j_{Tot}}$',interpreter='latex',location='best',fontsize=12);
leg1.ItemTokenSize(1)=10;
colormap(cmapAC)
axis tight
ylim([0 2.5])
axis square

t2=nexttile();
j_mean=plot(DOC, jTot_meanD2,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
    plot(DOC, jTot_allD2(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('DOC (\mug C l^{-1})')
ylabel('jTot')
ylim([0 2.5])

cbar=colorbar;
colormap(cmapAC)
ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square

t3=nexttile();
j_mean=plot(L, jTot_meanD3,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(L, jTot_allD3(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('L (\muE)')
ylabel('jTot')
ylim([0 2.5])

colormap(cmapAC)
ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square

t4=nexttile();
j_meanD=plot(Si, jTot_meanD4,'color',cmap(3,:),'lineStyle','-.' ,'linewidth',2);
hold on
for i=1:n
plot(Si, jTot_allD4(:,i),'color',cmapAC(i,:), 'linewidth',1)
end
xlabel('Si (\mug Si l^{-1})')
ylabel('jTot')
% colorbar
cbar=colorbar;
colormap(cmapAC)

ylabel(cbar,'Cell mass (\mug C)','FontSize',10)
cbar.Ticks=linspace(.1,1,n-1);%[.1 .5 1]; % correspond to size-classes 1,5,10
% pow1=ceil(log10(m(1)));
cbar.TickLabels={'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
axis square
axis tight
ylim([0 2.5])

sgtitle('Diatoms')


%% --------------------------------------------------------------------------
%--------------------------------------------------------------------------
function ProdNet = getProdNetGeneralists(this, u)
    % Precompute parts of the respiration term that don't depend on `i`
    resp = fTemp2 * this.Jresp + ...
           bL * this.JLreal + ...
           bN * this.JNreal .* this.JLreal ./ (this.JLreal + this.JDOCreal) + ...
           bg * this.Jnet .* this.JLreal ./ (this.JLreal + this.JDOCreal + this.JFreal);

    % Vectorized calculation of ProdNet
    ProdNet = sum(max(0, (this.JLreal - resp) .* u ./ this.m));
end
function ProdBact = getProdBactGeneralists(this, u)
    % Vectorized calculation of ProdBact
    ProdBact = sum(max(0, this.JDOC - fTemp2 * this.Jresp) .* u ./ this.m);
end
