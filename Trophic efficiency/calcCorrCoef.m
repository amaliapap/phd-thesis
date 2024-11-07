mte=squeeze(sim.ProdHTLAnnual./sim.ProdNetAnnual);
mte_vec=double(mte(:));
%%
T_arr=squeeze(mean(mean(sim.T,4,'omitnan'),1,'omitnan'));
L_arr=squeeze(mean(mean(sim.L,4,'omitnan'),1,'omitnan'));

T_vec=double(T_arr(:));
L_vec=double(L_arr(:));
TL_arr=[T_vec';L_vec']';
%%
[R,P]=corrcoef(TL_arr);