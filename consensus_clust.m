function [] = run_consensus_clust(hclust_dir, ccut)
%ccut = 0.7;
list_fns = dir([hclust_dir '\r_k*.mat']);
sz2=size(list_fns);
for j = 1:sz2(1)
    load([hclust_dir '\' list_fns(j).name],'r');
    size(r)
    r(r<ccut) = 0;
    S = eventSamples(r, 10000);
    [Sc, Tree] = hierarchicalConsensus(S);
    C = coclassificationMatrix(S);
    [ax_C,ax_H,o] = consensusPlot(C, Sc, Tree);
    [partitions,part_thr] = allPartitions(Sc,Tree);
    c=colorbar(ax_C)
    c.Location = 'westoutside'
    fn = [hclust_dir '\ccut' num2str(ccut) '_Coclass_' list_fns(j).name]
    save(fn, 'C')
    fn = [hclust_dir '\ccut' num2str(ccut) '_Sc_' list_fns(j).name]
    save(fn, 'Sc')
    fn = [hclust_dir '\ccut' num2str(ccut) '_Tree_' list_fns(j).name]
    save(fn, 'Tree')
    fn = [hclust_dir '\ccut' num2str(ccut) '_Partition_' list_fns(j).name]
    save(fn, 'partitions')
    fn = [hclust_dir '\ccut' num2str(ccut) '_Thr_Partition_' list_fns(j).name]
    save(fn, 'part_thr')
    fn = [hclust_dir '\ccut' num2str(ccut) '_o_Fig_' list_fns(j).name]
    save(fn, 'o')
    fn = [hclust_dir '\ccut' num2str(ccut) '_Fig_' list_fns(j).name '.png']
    saveas(ax_C,fn)
    fn = [hclust_dir '\ccut' num2str(ccut) '_Fig_' list_fns(j).name '.eps']
    saveas(ax_C,fn, 'epsc')
end
