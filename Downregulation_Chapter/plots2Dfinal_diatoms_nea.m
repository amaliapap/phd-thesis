% plots2D
% Varies the potential resource uptakes jL,jDOC,jN,jSi for values within
% the range of logNos
% and plots all the effective uptakes jX and downregulation factors dX
% for each varying jX
% The varying jXs are selected based on their influence on the monotonicity
% of the effective uptakes with respect to the varying potential uptake
% and the values that are kept constant are chosen based on the same
% criterion, using the function findRegimeShifts
% Calls findRegimeShifts.m
clc; clear all;
%%
% colormaps
alph3=0.6*ones(3,1);
alph4=0.6*ones(4,1);

threecolores=brewermap(3,'Dark2');
fourcolores=brewermap(4,'Dark2');
cmap3=cat(2,threecolores,alph3);
cmap4=cat(2,fourcolores,alph4);

%%
this.n=10;
p=parametersDiatoms(this.n);

Jmax = p.alphaJ * p.m .* (1.0 - p.nu)./p.m; % max rate (1/day)
sizeClass=4;
m=p.m(sizeClass);
%
% isPureAutotroph=false; % or false to include phagotrophy
withAffinities=false;

logNos=logspace(-3,0,10);
% logNos=linspace(0.01,1,10);

jDOC=logNos;
jN=logNos;
jL=logNos;
jSi=logNos;
JXgradients=[jL; jDOC; jN; jSi]';

var_i=jL;
var_i_str='j_L';
var_j=jDOC;
var_j_str='j_{DOC}';
var_k=jN;
var_k_str='j_N';
var_l=jSi;
var_l_str='j_{Si}';
var_str={var_i_str; var_j_str; var_k_str; var_l_str};
% var_rank shows the most influential jX (when varied) in descending order (index based on their position in JXgradients)
% table2(i,j): 
[var_rank,table2]=findRegimeShifts(true,JXgradients);
most_frequent_index=cell2mat(mostFrequentElement(table2)); % returns the most frequent element in each row
influenc1=var_str(var_rank(1));
influenc2=var_str(var_rank(2));
%% Changing jL
for i=1:length(var_i)
    for j=5%most_frequent_index(2)%1:length(var_j)
        for k=5%most_frequent_index(3)%1:length(var_k)
            for l=6%most_frequent_index(4)%1:length(var_l)
                thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities);
                thisL{i,j,k,l}=costDiatoms(thisRate);
                JLmon(i)=thisL{i,j,k,l}.JLreal;
                JNmon(i)=thisL{i,j,k,l}.JN;
                JDOCmon(i)=thisL{i,j,k,l}.JDOC;
                JSimon(i)=thisL{i,j,k,l}.JSi;
                Jtot(i)=thisL{i,j,k,l}.Jtot;
                dL(i)=thisL{i,j,k,l}.dL;
                dN(i)=thisL{i,j,k,l}.dN;
                dDOC(i)=thisL{i,j,k,l}.dDOC;
                dSi(i)=thisL{i,j,k,l}.dSi;
            end
        end
    end
end
%  Varying jN
changing_var=var_k_str;
for i=7%most_frequent_index(1)%1:length(var_i)
    for j=5%most_frequent_index(2)%1:length(var_j)
        for k=1:length(var_k)
            for l=6%most_frequent_index(4)%1:length(var_l)
                thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities);
                thisL{i,j,k,l}=costDiatoms(thisRate);
                JLmon2(k)=thisL{i,j,k,l}.JLreal;
                JNmon2(k)=thisL{i,j,k,l}.JN;
                JDOCmon2(k)=thisL{i,j,k,l}.JDOC;
                JSimon2(k)=thisL{i,j,k,l}.JSi;
                Jtot2(k)=thisL{i,j,k,l}.Jtot;
                dL2(k)=thisL{i,j,k,l}.dL;
                dN2(k)=thisL{i,j,k,l}.dN;
                dDOC2(k)=thisL{i,j,k,l}.dDOC;
                dSi2(k)=thisL{i,j,k,l}.dSi;
            end
        end
    end
end
%% PLOTS
x0=1;
y0=2;
height=9;  
width=16;

fig_no=1;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Varying fluxes Diatims: jL,jN');

clf
set(gcf,'color','w');
tiledlayout(2,2,Padding="compact",TileSpacing="tight",TileIndexing="columnmajor") 
set(groot,'defaultlinelinewidth',2)
 plot2tilesJX(var_i_str,jL,JLmon,JDOCmon,JNmon,JSimon,Jtot,dL,dDOC,dN,dSi)
 plot2tilesJX(var_k_str,jN,JLmon2,JDOCmon2,JNmon2,JSimon2,Jtot2,dL2,dDOC2,dN2,dSi2,true,true);
 
%% Changing jL
for i=1:length(var_i)
    for j=5%most_frequent_index(2)%1:length(var_j)
        for k=5%most_frequent_index(3)%1:length(var_k)
            for l=5%most_frequent_index(4)%1:length(var_l)
                thisRate=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities);
                thisL{i,j,k,l}=costGeneralists(thisRate);
                JLmon(i)=thisL{i,j,k,l}.JLreal;
                JNmon(i)=thisL{i,j,k,l}.JN;
                JDOCmon(i)=thisL{i,j,k,l}.JDOC;
                JSimon(i)=thisL{i,j,k,l}.JFreal;
                Jtot(i)=thisL{i,j,k,l}.Jtot;
                dL(i)=thisL{i,j,k,l}.dL;
                dN(i)=thisL{i,j,k,l}.dN;
                dDOC(i)=thisL{i,j,k,l}.dDOC;
            end
        end
    end
end
%  Varying jN
changing_var=var_k_str;
for i=7%most_frequent_index(1)%1:length(var_i)
    for j=5%most_frequent_index(2)%1:length(var_j)
        for k=1:length(var_k)
            for l=5%most_frequent_index(4)%1:length(var_l)
                thisRate=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities);
                thisL{i,j,k,l}=costGeneralists(thisRate);
                JLmon2(k)=thisL{i,j,k,l}.JLreal;
                JNmon2(k)=thisL{i,j,k,l}.JN;
                JDOCmon2(k)=thisL{i,j,k,l}.JDOC;
                JSimon2(k)=thisL{i,j,k,l}.JFreal;
                Jtot2(k)=thisL{i,j,k,l}.Jtot;
                dL2(k)=thisL{i,j,k,l}.dL;
                dN2(k)=thisL{i,j,k,l}.dN;
                dDOC2(k)=thisL{i,j,k,l}.dDOC;
            end
        end
    end
end
%% PLOTS
x0=1;
y0=2;
height=9;  
width=16;

fig_no=2;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Varying fluxes Generalists: jL,jN');

clf
set(gcf,'color','w');
tiledlayout(2,2,Padding="compact",TileSpacing="tight",TileIndexing="columnmajor") 
set(groot,'defaultlinelinewidth',2)
 plot2tilesJX(var_i_str,jL,JLmon,JDOCmon,JNmon,JSimon,Jtot,dL,dDOC,dN,0,false)
 plot2tilesJX(var_k_str,jN,JLmon2,JDOCmon2,JNmon2,JSimon2,Jtot2,dL2,dDOC2,dN2,0,false,true)

%%
% ax = gcf;
% exportgraphics(ax,'diatomsVarJLjN.pdf')

function plot2tilesJX(var_i_str,xval,yval1,yval2,yval3,yval4,yval5,yval6,yval7,yval8,yval9,isDiatom,rightPlot)
arguments
    var_i_str 
    xval 
    yval1 
    yval2 
    yval3 
    yval4 
    yval5 
    yval6 
    yval7 
    yval8 
    yval9 = 0;
    isDiatom  = true
    rightPlot = false;
end
fourcolores=brewermap(4,'Dark2');
threecolores=brewermap(3,'Dark2');

t1=nexttile; % Varying jL
[~,idx_a]=find(yval1==yval2); % find regime-shift points
[~,idx_b]=find(yval3==yval2);

for i=1:length(xval)
    plot(xval,yval1)
    hold on
    plot(xval,yval2,'o')
    plot(xval,yval3)
    plot(xval,yval4,':')
    plot(xval,yval5,'k:',LineWidth=2.5)
end
if ~isempty(idx_a)
xline(xval(min(idx_a)),'k--')
end
if ~isempty(idx_b)
xline(xval(min(idx_b)),'k--')
end

% axis square
if rightPlot==false
    leg2=legend('j_{Lreal}','j_{DOC}','j_N','j_{Si}','j_{tot}');
    if isDiatom==false
        leg2=legend('j_{Lreal}','j_{DOC}','j_N','j_{F}','j_{tot}');
    end
    % xlabel(append(var_i_str,' (day^{-1})'))
    ylabel('Rates (day^{-1})')
    leg2.ItemTokenSize(1)=10;
    leg2.Location='northwestoutside';
    plotlabel('a',false)
else
    plotlabel('b',false)
end
    colororder(t1,[fourcolores; 0,0,0]);
    box on
set(gca,'XScale','log')
axis tight
%
% delta's
t2=nexttile;

for i=1:length(xval)
    hold on
    plot(xval,yval6)
    plot(xval,yval7,'o')
    plot(xval,yval8)
    if isDiatom==true
        plot(xval,yval9,':')
    end
end
if ~isempty(idx_a)
    xline(xval(min(idx_a)),'k--')
end
if ~isempty(idx_b)
    xline(xval(min(idx_b)),'k--')
end
% axis square
if isDiatom==false
    colororder(t2,threecolores);
else
    colororder(t2,fourcolores);
end
if rightPlot==false
    if isDiatom==false
        leg1=legend('\delta_L','\delta_{DOC}','\delta_N');
        % colororder(t2,threecolores);
    else
        leg1=legend('\delta_L','\delta_{DOC}','\delta_N','\delta_{Si}');
        % colororder(t2,fourcolores);
    end
    ylabel('deltas (-)')
    leg1.ItemTokenSize(1)=10;
    leg1.Location='northwestoutside';
    plotlabel('c',false)
else
    plotlabel('d',false)
end
xlabel(append(var_i_str,' (day^{-1})'))
    box on

set(gca,'XScale','log')
end