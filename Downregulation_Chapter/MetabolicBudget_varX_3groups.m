% Calls: calcThisDiatom and costDiatoms
%
%           [jLcost;    jNcost;   jSicost;   jDOCcost; jNet_cost; jR_mean;  Jnet;      Jpassive]
% newcolors0 = {'#b5fded',' #b5eded','#b5dded','#b5cded','#b5bded', '#d7bde2','#b25781', '#66b9bb'};
% newcolors = {'#b5fded','#8ef2ff','#fffad3','#c5b1e4','#b5bded', '#d7bde2','#46c6c6','#d7a1b8'};
% seasonality
cmap=flip(cmocean('deep',7));
cmap3=cmap(2:end,:);
cmap4=(cmocean('matter',6));
ccmap2=[cmap3(1:end,:);cmap4(4:end-1,:)];

% blue-purple
% cmap=brewermap(21,'BrBG');
% cmap3=cmap(13:18,:);
% cmap4=brewermap(10,'PRGn');
% ccmap=cmap4(3:6,:);
% ccmap2=[cmap3(1:end-2,:);[215, 219, 221]/250;ccmap];
newcolors=flip(ccmap2);
% newcolors = [ccmap2(3,:);ccmap2(4,:);ccmap2(5,:);ccmap2(2,:);ccmap2(7,:);ccmap2(1,:);[215, 219, 221]/250; [104, 186, 187]/250];
% Calculate Jmax

this.n=10;
p=parametersDiatoms(this.n);

Jmax = p.alphaJ * p.m .* (1.0 - p.nu)./p.m; % max rate (1/day)
sizeClass=4;
m=p.m(sizeClass);
%% Diatom downregulation and synthesis
withAffinities=false;
specificRates=true;

logNos=linspace(0.1,1,5);
% logNos=logspace(-3,0,5);


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
%% Changing jL
for i=1:length(var_i)
    for j=2%1:length(var_j)
        for k=2%1:length(var_k)
            for l=2%1:length(var_l)
            thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities,specificRates);
            thisL=costDiatoms(thisRate);
            thisRateg=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,false,withAffinities,specificRates);
            thisG=costGeneralists(thisRateg);
            thisRatep=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,true,withAffinities,specificRates);
            thisP=costGeneralists(thisRatep);
            
            varL.jResptot(:,i)=thisL.Jresptot;
            varL.MB(:,i)=thisL.MB;
            varL.MBdetails(:,i)=thisL.MBdetails;

            varLg.jResptot(:,i)=thisG.Jresptot;
            varLg.MB(:,i)=thisG.MB;
            varLg.MBdetails(:,i)=thisG.MBdetails;
            varLp.jResptot(:,i)=thisP.Jresptot;
            varLp.MB(:,i)=thisP.MB;
            varLp.MBdetails(:,i)=thisP.MBdetails;
            end
        end
    end
end

%% Changing jN
for i=3%length(var_i)
    for j=2%1:length(var_j)
        for k=1:length(var_k)
            for l=2%1:length(var_l)
            thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities,specificRates);
            thisL=costDiatoms(thisRate);
            thisRateg=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,false,withAffinities,specificRates);
            thisG=costGeneralists(thisRateg);
            thisRatep=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,true,withAffinities,specificRates);
            thisP=costGeneralists(thisRatep);
            
            varN.jResptot(:,k)=thisL.Jresptot;
            varN.MB(:,k)=thisL.MB;
            varN.MBdetails(:,k)=thisL.MBdetails;
            varNg.jResptot(:,k)=thisG.Jresptot;
            varNg.MB(:,k)=thisG.MB;
            varNg.MBdetails(:,k)=thisG.MBdetails;
            varNp.jResptot(:,k)=thisP.Jresptot;
            varNp.MB(:,k)=thisP.MB;
            varNp.MBdetails(:,k)=thisP.MBdetails;
            end
        end
    end
end
%% Changing jDOC
for i=3%length(var_i)
    for j=1:length(var_j)
        for k=2%1:length(var_k)
            for l=2%1:length(var_l)
            thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities,specificRates);
            thisL=costDiatoms(thisRate);
            thisRateg=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,false,withAffinities,specificRates);
            thisG=costGeneralists(thisRateg);
            thisRatep=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,true,withAffinities,specificRates);
            thisP=costGeneralists(thisRatep);
            
            varDOC.jResptot(:,j)=thisL.Jresptot;
            varDOC.MB(:,j)=thisL.MB;
            varDOC.MBdetails(:,j)=thisL.MBdetails;
            varDOCg.jResptot(:,j)=thisG.Jresptot;
            varDOCg.MB(:,j)=thisG.MB;
            varDOCg.MBdetails(:,j)=thisG.MBdetails;
            varDOCp.jResptot(:,j)=thisP.Jresptot;
            varDOCp.MB(:,j)=thisP.MB;
            varDOCp.MBdetails(:,j)=thisP.MBdetails;
            end
        end
    end
end

%% Changing jSi and jF
for i=3%length(var_i)
    for j=2%1:length(var_j)
        for k=2%1:length(var_k)
            for l=1:length(var_l)
            thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities,specificRates);
            thisL=costDiatoms(thisRate);
            thisRateg=calcThisGeneralistX(jN(k),jL(i),jSi(l),jDOC(j),m,false,withAffinities,specificRates);
            thisG=costGeneralists(thisRateg);
           
            varSi.jResptot(:,l)=thisL.Jresptot;
            varSi.MB(:,l)=thisL.MB;
            varSi.MBdetails(:,l)=thisL.MBdetails;
            varFg.jResptot(:,l)=thisG.Jresptot;
            varFg.MB(:,l)=thisG.MB;
            varFg.MBdetails(:,l)=thisG.MBdetails;
            
            end
        end
    end
end

%%
x0=1;
y0=2;
height=14;  
width=13;

fig_no=11;
fig=figure(fig_no);
set(fig,'Renderer','Painters','Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto','Name','Metabolic budget cell');

clf
set(gcf,'color','w');
tiledlayout(4,3,Padding="tight",TileIndexing="columnmajor",TileSpacing="tight") % column1: eutrophic  column2: oligotrophic
%
plotRates=false;
% Phototrophs
t1=nexttile(1);
    plotNexttiles(var_i,varLp,'j_{L} (day^{-1})',plotRates)
    title(t1,'Phototrophs')
plotlabel('a',false)
t2=nexttile(2);
    plotNexttiles(var_k,varNp,'j_{N} (day^{-1})',plotRates)
    plotlabel('d',false)

t3=nexttile(3);
    plotNexttiles(var_j,varDOCp,'j_{DOC} (day^{-1})',plotRates)
    plotlabel('g',false)

%
% Diatoms
td1=nexttile(1+4);
    plotNexttiles(var_i,varL,'j_{L} (day^{-1})',plotRates)
    title(td1,'Diatoms')
    plotlabel('b',false)
td2=nexttile(2+4);
    plotNexttiles(var_k,varN,'j_{N} (day^{-1})',plotRates)
    plotlabel('e',false)
td3=nexttile(3+4);
    plotNexttiles(var_j,varDOC,'j_{DOC} (day^{-1})',plotRates)
    plotlabel('h',false)
td4=nexttile(4+4);
    plotNexttiles(var_l,varSi,'j_{Si} (day^{-1})',plotRates)
    plotlabel('j',false)
%
% Generalists
tg1=nexttile(1+8);
    plotNexttiles(var_i,varLg,'j_{L} (day^{-1})',plotRates)
    title(tg1,'Generalists')
    ylabel('(day^{-1})')
    plotlabel('c',false)
tg2=nexttile(2+8);
    plotNexttiles(var_k,varNg,'j_{N} (day^{-1})',plotRates)
    ylabel('(day^{-1})')
    plotlabel('f',false)
tg3=nexttile(3+8);
    plotNexttiles(var_j,varDOCg,'j_{DOC} (day^{-1})',plotRates)
    ylabel('(day^{-1})')
    plotlabel('i',false)
tg4=nexttile(4+8);
    plotNexttiles(var_l,varFg,'j_{F} (day^{-1})',plotRates)
    ylabel('(day^{-1})')
    plotlabel('k',false)
    leg3=legend('Basal','N cost','Si/F cost','DOC cost',...
    'Synth. cost','L cost','Jnet','Passive loss','Division rate','\SigmaJresp','Location','bestoutside');
    leg3.ItemTokenSize(1)=10;
    leg3.NumColumns=2;
    leg3.Position(1)=0.01;
    leg3.Position(2)=.08;
colororder(newcolors)
axx= findall(gcf, 'Type', 'Axes');
for i = 2:length(axx)
    set(axx(i), 'YLim', [0 0.7]);
    axx(i).YLim  = [0 0.7];
end


%%
function plotNexttiles(var_i,varL,str,plotRates)
arguments
    var_i 
    varL 
    str 
    plotRates=false
end
cmap=brewermap(21,'BrBG');
cmap3=cmap(13:18,:);
% palePink=cmap3(2,:);
cmap4=brewermap(10,'PRGn');
ccmap=cmap4(3:6,:);
ccmap2=[cmap3(1:end-2,:);[215, 219, 221]/250;ccmap];
% newcolors=flip(ccmap2);

if plotRates==false
    % bar(var_i,varL.MBdetails./sum(varL.MBdetails,1),'stacked','FaceAlpha',0.7)
    panelContinuousBars(var_i,varL,false)
    yyaxis right
    line_col=[215, 219,221]/250;
    plot(var_i,varL.MB(1,:),'Color',line_col,LineStyle='-.',LineWidth=2)
    hold on
    plot(var_i,varL.jResptot,'Color',line_col,LineWidth=2)

    ax=gca;
    ax.YAxis(1).Color = '#b25782';
    ax.YAxis(2).Color = [.4 .4 .4];

    ax.XTick=var_i;
    xtickformat('%.1f')

else
    plot(var_i,varL.MBdetails)
    ylim([0, 1])
end
% colororder(newcolors);
xlabel(str)
axis tight
end