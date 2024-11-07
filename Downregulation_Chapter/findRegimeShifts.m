% Finds values for the fluxes that induce regime shifts when one variable
% changes at each time- these values will then be used for
% calls: calcRegimeShifts.m
% Input: JXgradients=[jL; jDOC; jN; jSi]'  is a N-by-4 matrix 
% 
% Outputs:
% - var_rank: shows the most influential jX (when varied) in descending
%             order (index based on their position in JXgradients)
% - table2: 
% Example:
% table2 =
%    JLmon JDOCmon JNmon JSimon (=dy)
%jL    7       1     1     3
%jDOC  2       1     1     1
%jN    1       1     1     1
%jSi   1       1     1     2
%(dx)
% e.g. to achieve the biggest change in monotonicity for JLmon
% the potential uptakes should have the values corresponding 
% to their index in JXgradients. jL =JXgradients(7,1)
function [z,tableJXbyJXmon]=findRegimeShifts(isDiatom,JXgradients)
%%
% We need to see tableJXbyJXmon and decide on which indices to use
% isDiatom=true;
%%
this.n=10;
if isDiatom==true
    p=parametersDiatoms(this.n);
else
    p=parametersGeneralists(this.n);
end
Jmax = p.alphaJ * p.m .* (1.0 - p.nu)./p.m;
sizeClass=4;
m=p.m(sizeClass);
%%
isPureAutotroph=false; % or false to include phagotrophy
withAffinities=false;

% %% Changing all variables
% jVar = linspace(0,Jmax(sizeClass),10);
% jDOC=jVar;
% jN=jVar;
% jL=jVar;
% jSi=jVar;
% 
% JXgradients=[jL; jDOC; jN; jSi]';
jL=JXgradients(:,1);
jDOC=JXgradients(:,2);
jN=JXgradients(:,3);
jSi=JXgradients(:,4);

var_i=jL;
var_i_str='j_L';
var_j=jDOC;
var_j_str='j_{DOC}';
var_k=jN;
var_k_str='j_N';
var_l=jSi;
var_l_str='j_{Si}';
% changing_var=var_k_str;
for i=1:length(var_i)
    for j=1:length(var_j)
        for k=1:length(var_k)
            for l=1:length(var_l)
                if isDiatom==true
                    thisRate=calcThisDiatomsX(jN(k),jL(i),jSi(l),jDOC(j),m,withAffinities);
                    thisL{i,j,k,l}=costDiatoms(thisRate);
                    JSimon(i,j,k,l)=thisL{i,j,k,l}.JSi;
                else
                    jF=jSi;
                    thisRate=calcThisGeneralistX(jN(k),jL(i),jF(l),jDOC(j),m,isPureAutotroph,withAffinities);
                    thisL{i,j,k,l}=costGeneralists(thisRate);
                    JSimon(i,j,k,l)=thisL{i,j,k,l}.JF;
                end
                JNmon(i,j,k,l)= thisL{i,j,k,l}.JN;
                JDOCmon(i,j,k,l)=thisL{i,j,k,l}.JDOC;
                JLmon(i,j,k,l)=thisL{i,j,k,l}.JLreal;
            end
        end
    end
end

%%
% Allocate all effective uptakes J'X'mon in JXmon:(4x10x10x10x10)
nv=length(JXgradients);
JXmon=zeros(4,nv,nv,nv,nv);
JXmon(1,:,:,:,:)=JLmon;
JXmon(2,:,:,:,:)=JDOCmon;
JXmon(3,:,:,:,:)=JNmon;
JXmon(4,:,:,:,:)=JSimon;

idxJXmonJX=zeros(4,4,4);
max_dJXmon_dJX=zeros(4);

for h=1:4 %JXmon
    JX_mon=squeeze(JXmon(h,:,:,:,:));
    for i=1:4 %
        [maxVal,idx]=calcRegimeShifts(JX_mon,i,JXgradients);
        idxJXmonJX(:,h,i)=idx;
        max_dJXmon_dJX(h,i)=maxVal'; % maximum change in monotonicity for denominator JXgardients(:,i)
    end
end
% find the column of the most influential jX for all JXmon on average
JXmonJX_mean=mean((max_dJXmon_dJX));
[~,jX_i]=max(mean(max_dJXmon_dJX,1)); 

x=JXmonJX_mean;
y=sort(x,'descend'); 
[~,Locb] = ismember(y,x);
z=Locb'; % returns the most influential values (most to least)

% Array of optimal indices corresponding to max monotonicity-change ...
% (having the most influential jX as denominator);each column corresponds to JXmon 
tableJXbyJXmon=squeeze(idxJXmonJX(:,:,jX_i));

%% Same as ln:72-79 without changing JXmon
% idxJXmonJX=zeros(4);
% 
% JXgradients=[jL; jDOC; jN; jSi]';
% JX_mon=squeeze(JXmon(4,:,:,:,:));
% for i=1:4 %
%     [maxVal,idx]=calcRegimeShifts(JX_mon,i,JXgradients);
%     idxJXmonJX(:,i)=idx;
%     max_dJXmon_dJX(i)=maxVal'; % maximum change in monotonicity for denominator JXgardients(:,i)
% end
%
% explantation: max(maxJLmonJX)=4, means that JL has the biggest change in
%               monotoicity as a function of the 4th gradient which is jSi
% the other fluxes should be idxJLmonJX(1:3,4)


% %% This is the algorithm used in calcRegimeShifts.m
% for i=1:length(var_i)
%     for j=1:length(var_j)
%         for k=1:length(var_k)
%             for l=1:length(var_l)  %        This calculates the gradient with respect to resource N
%                 dydxJL4(i,j,:,l)= gradient(squeeze(JLmon(i,j,:,l))) ./ gradient(JXgradients(1:10,3));
%                 gradient_mean4(i,j,k,l)=mean(abs(dydxJL4(i,j,k,l)));
%             end
%         end
%     end
% end
% [mxv4,idx4] = max(gradient_mean4(:));
% [r,c,w,p] = ind2sub(size(gradient_mean),idx4);
% idx_maxJN=[r,c,w,p];

% Investigates the monotonicity of the curves: (JXmon/JXgradients(:,dim_JXgradients)
%
% example: if dim_JXgradients=3 and JXmon=JLmon then it returns the indices
%          that produce the sharper changes in monotoncity of JLmon 
%          with respect to JN(prescribed) 
% 
function [mxv4,idx_max]=calcRegimeShifts(JXmon,dim_JXgradients,JXgradients)

var_i=JXgradients(:,1);
% var_i_str='j_L';
var_j=JXgradients(:,2);
% var_j_str='j_{DOC}';
var_k=JXgradients(:,3);
% var_k_str='j_N';
var_l=JXgradients(:,4);
% var_l_str='j_{Si}';
% gradient_mean4=zeros(length(var_i),length(var_j),length(var_k),length(var_l));
% dydxJXmonJx=gradient_mean4;
%% 
% dim_JXgradients=4;
% JXmon=JLmon;
for i=1:length(var_i)
    ii=i;
    for j=1:length(var_j)
        jj=j;
        for k=1:length(var_k)
            kk=k;
            for l=1:length(var_l)
                ll=l;
                if dim_JXgradients==1
                    ii=1:length(var_i);
                elseif dim_JXgradients==2
                    jj=1:length(var_j);
                elseif dim_JXgradients==3
                    kk=1:length(var_k);
                elseif dim_JXgradients==4
                    ll=1:length(var_l);
                end
                % This calculates the gradient of JL mon with respect to resource N
                tempval= gradient(squeeze(JXmon(ii,jj,kk,ll)))...
                    ./ gradient(JXgradients(:,dim_JXgradients));
                if size(tempval,2)>1
                    tempval= gradient(squeeze(JXmon(ii,jj,kk,ll)))...
                    ./ gradient(JXgradients(:,dim_JXgradients)');
                end
                % size1=size(gradient(squeeze(JXmon(ii,jj,kk,ll))))
                % size2=size(gradient(JXgradients(:,dim_JXgradients)))
                % size(tempval)
                dydxJXmonJx(ii,jj,kk,ll)=tempval;
                gradient_mean4(i,j,k,l)=mean(abs(dydxJXmonJx(i,j,k,l)));
            end
        end
    end
end
[mxv4,idx4] = max(gradient_mean4(:));
[r,c,w,p] = ind2sub(size(gradient_mean4),idx4);
idx_max=[r,c,w,p];


