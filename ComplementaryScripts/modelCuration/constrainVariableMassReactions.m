%
% FILE NAME:    constrainVariableMassReactions.m
% 
% DATE CREATED: 2018-09-28
%     MODIFIED: 2018-10-02
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to identify and constrain all reactions that involve the
%          same reactants producing different compounds, such that the
%          different compounds vary in mass. For example:
%
%           LCAT39e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[s]
%            LCAT5e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[s]
%           LCAT31e: cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[s]
%
%          These reactions have the same reactants (cholesterol and PC-LD
%          pool), but differ in only one product, which results in a PC-LD
%          pool of variable mass (since some of these reactions are
%          reversible, and some of the mets involved can be broken down
%          into their components).


%% Script for identifying all potential reaction sets

% if ~exist('ihuman','var')
%     load('humanGEM.mat');
% end
% 
% massVariable_rxn_sets = [];
% massConsistent_rxn_sets = [];
% 
% for n = 4:100
%     
%     % find all reactions with X metabolites
%     rxn_ind = sum(ihuman.S ~= 0,1)' == n;
%     
%     % find all mets involved in these reactions
%     met_ind = any(ihuman.S(:,rxn_ind) ~= 0,2);
%     
%     % extract subset of stoich matrix for these reactions and metabolites
%     s = ihuman.S(met_ind,rxn_ind);
%     
%     % keep track of original reaction indices
%     orig_rxn_ind = (1:length(ihuman.rxns))';
%     orig_rxn_ind = orig_rxn_ind(rxn_ind);
%     
%     % calculated the hamming distance between these reactions
%     rxn_dist = squareform(pdist(s','hamming'));
%     
%     % we are interested in reactions that differ by only one metabolite, which
%     % would correspond to a hamming distance of 2/Nmets
%     d = 2/sum(met_ind);
%     
%     % find all cases where at least 2 other rxns differ by this distance
%     check_rxns = find(sum(rxn_dist == d) >= 2);
%     
%     % obtain unique list of similar reaction sets (do this by retrieving rxns
%     % with the specified hamming distance, as well as a distance of zero, so
%     % that the reaction set includes the checked reaction itself).
%     rxn_sets = unique(ismember(rxn_dist(check_rxns,:),[0,d]),'rows');
%     
%     
%     % remove any sets that are just a subset of another
%     i = 0;
%     while ~isempty(rxn_sets)
%         i = i+1;
%         r = rxn_sets - rxn_sets(i,:);
%         if sum(all(r >= 0,2)) > 1
%             rxn_sets(i,:) = [];
%             i = 0;
%         elseif i == size(rxn_sets,1)
%             break
%         end
%     end
%         
%     % convert to logical that corresponds to original rxn indexing
%     rxn_sets = rxn_sets';  % transpose matrix
%     rxn_sets_orig = false(length(ihuman.rxns),size(rxn_sets,2));
%     for i = 1:size(rxn_sets,2)
%         rxn_sets_orig(orig_rxn_ind(rxn_sets(:,i)),i) = true;
%     end
%     
%     % for each reaction set, determine which differ by metabolites that vary in
%     % their formula
%     ignore_sets = [];
%     for i = 1:size(rxn_sets_orig,2)
%         diff_mets = sum(ihuman.S(:,rxn_sets_orig(:,i)) ~= 0, 2) == 1;
%         if length(unique(ihuman.metFormulas(diff_mets))) == 1
%             ignore_sets = [ignore_sets; i];
%         end
%     end
%     
%     % append these reaction sets to the existing sets
%     if ~isempty(ignore_sets)
%         massConsistent_rxn_sets = [massConsistent_rxn_sets, rxn_sets_orig(:,ignore_sets)];
%         rxn_sets_orig(:,ignore_sets) = [];
%     end
%     massVariable_rxn_sets = [massVariable_rxn_sets, rxn_sets_orig];
%     
% end
% 
% massVariable_rxn_sets = logical(massVariable_rxn_sets);
% massConsistent_rxn_sets = logical(massConsistent_rxn_sets);


%% Constrain reaction sets after manual curation of results from above
% After manually curating the reaction sets in "massVariable_rxn_sets" to
% confirm those that violate mass balances, the following list of reactions
% will be constrained to zero, to avoid creation or destruction of mass or
% energy:
%
% RXN ID    RXN EQUATION
%
% AGPAT1    2 H+[c] + R Total 2 Coenzyme A[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT2    palmitoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT3    oleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
% AGPAT4    linoleoyl-CoA[c] + Lysophosphatidic Acid[c] => CoA[c] + phosphatidate-LD-TAG pool[c]
%  
% LCAT55e   cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphoinositol[s]
% LCAT17e   cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphoinositol[s]
% LCAT4e    cholesterol[s] + PI pool[s] => cholesterol-ester pool[s] + 1-Arachidonoylglycerophosphoinositol[s]
%  
% LCAT12e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Oleoylglycerophosphoethanolamine (Delta 9)[s]
% LCAT54e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphoethanolamine[s]
% LCAT16e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphoethanolamine[s]
% LCAT19e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 2-Linoleoylglycerophosphoethanolamine[s]
% LCAT3e    cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Arachidonoyl-Sn-Glycero-3-Phosphoethanolamine[s]
% LCAT40e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatrienoylglycerophosphoethanolamine (Delta 11, 14, 17), Lpe (20:3)[s]
% LCAT41e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglyceroethanolamine (Delta 4, 7, 10, 13, 16, 19), Lpe (22:6)[s]
% LCAT42e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Docosatetraenoyglycerophosphoethanolamine (22:4, Delta 7, 10, 13, 16)[s]
% LCAT43e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Dihomo-Linolenoylglycerophosphoethanolamine (20:3, Delta 8, 11, 14)[s]
% LCAT44e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Didecanoylglycerophosphoethanolamine (C12:0 Pe)[s]
% LCAT45e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Myristoylglycerophosphoethanolamine (C14:0 Pe)[s]
% LCAT56e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Hexadecenoylglycerophosphoethanolamine (C16:1 Pe, Delta 9)[s]
% LCAT46e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Tridecanoylglycerophosphoethanolamine (C13:0 Pe)[s]
% LCAT47e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Pentadecanoylglycerophosphoethanolamine (C15:0 Pe)[s]
% LCAT48e   cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Heptadecanoylglycerophosphoethanolamine (C17:0 Pe)[s]
% LCAT9e    cholesterol[s] + PE-LD pool[s] => cholesterol-ester pool[s] + 1-Linoleoylglycerophosphoethanolamine (Delta 9,12)[s]
% 
% LCAT10e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Myristoylglycerophosphocholine[s]
% LCAT11e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Oleoylglycerophosphocholine (Delta 9)[s]
% LCAT13e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoleoylglycerophosphocholine (Delta 9)[s]
% LCAT14e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Palmitoylglycerophosphocholine[s]
% LCAT15e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Stearoylglycerophosphocholine[s]
% LCAT18e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Linoleoylglycerophosphocholine[s]
% LCAT20e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Oleoylglycerophosphocholine[s]
% LCAT21e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Palmitoylglycerophosphocholine[s]
% LCAT22e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 2-Stearoylglycerophosphocholine[s]
% LCAT23e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Gamma-Linolenoyl-Cholesterol, Cholesterol-Ester (18:3, Delta 6, 9, 12)[s]
% LCAT25e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Vaccenoyl-Cholesterol, Cholesterol-Ester (18:1, Delta 11)[s]
% LCAT26e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + 1-Timnodnoyl-Cholesterol, Cholesterol-Ester (20:5, Delta 5,8,11,14,17)[s]
% LCAT27e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + Cholesteryl Arachidonate, Cholesterol-Ester (20:4, Delta 5,8,11,14)[s]
% LCAT28e   cholesterol[s] + PC-LD pool[s] => 2-lysolecithin pool[s] + Cholesteryl Docosahexanoate, Cholesterol-Ester (22:6, Delta 4,7,10,13,16,19)[s]
% LCAT29e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Pentadecanoylglycerophosphocholine, Sn1-Lpc (15:0)[s]
% LCAT2e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Arachidonoyl-Glycero-3-Phosphocholine[s]
% LCAT30e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 9, 12, 15)[s]
% LCAT31e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[s]
% LCAT32e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Nonadecanoylglycerophosphocholine, Sn1-Lpc (19:0)[s]
% LCAT33e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosenoylglycerophosphocholine (Delta 11) ,Sn1-Lpc (20:1)[s]
% LCAT34e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatetraenoylglycerophosphocholine (Delta 8, 11, 14, 17), Sn1-Lpc (20:4)[s]
% LCAT35e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosapentenoylglycerophosphocholine (Delta 5, 8, 11, 14, 17), Sn1-Lpc (20:5)[s]
% LCAT36e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosatetraenoylglycerophosphocholine (Delta 7, 10, 13, 16), Sn1-Lpc (22:4)[s]
% LCAT37e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosapentenoylglycerophosphocholine (Delta 7, 10, 13, 16, 19), Sn1-Lpc (22:5)-W3[s]
% LCAT38e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosapentenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16), Sn1-Lpc (22:5)-W6[s]
% LCAT39e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[s]
% LCAT49e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Dihomo-Linolenoylglycerophosphocholine (20:3, Delta 8, 11, 14), Lysopc A C20:3[s]
% LCAT50e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Lignocericylglycerophosphocholine (24:0), Lysopc A C24[s]
% LCAT51e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C26:1 (Delta 5)[s]
% LCAT52e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C28:1 (Delta 5)[s]
% LCAT53e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + Lysopc A C28:0[s]
% LCAT57e   cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Docosahexaenoylglycerophosphocholine[s]
% LCAT5e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[s]
% LCAT6e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Eicosatrienoylglycerophosphocholine (Delta 11, 14, 17)[s]
% LCAT7e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Heptadecanoylglycerophosphocholine[s]
% LCAT8e    cholesterol[s] + PC-LD pool[s] => cholesterol-ester pool[s] + 1-Linoleoylglycerophosphocholine (Delta 9,12)[s]
% 
% SMS1      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/14:0), Sphingomyelin[c]
% SMS10     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/21:0), Sphingomyelin[c]
% SMS11     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/22:1), Sphingomyelin[c]
% SMS12     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/22:0), Sphingomyelin[c]
% SMS16     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/23:0), Sphingomyelin[c]
% SMS13     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/24:1), Sphingomyelin[c]
% SMS14     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/24:0), Sphingomyelin[c]
% SMS15     ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:0/25:0), Sphingomyelin[c]
% SMS2      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/15:0), Sphingomyelin[c]
% SMS3      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/16:1), Sphingomyelin[c]
% SMS4      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/16:0), Sphingomyelin[c]
% SMS5      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/17:0), Sphingomyelin[c]
% SMS6      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/18:0), Sphingomyelin[c]
% SMS7      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/18:1), Sphingomyelin[c]
% SMS8      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/20:1), Sphingomyelin[c]
% SMS9      ceramide pool[c] + PC-LD pool[c] => 1,2-diacylglycerol-LD-TAG pool[c] + Sm (D18:1/20:0), Sphingomyelin[c]
% 
% PEOLE_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Oleoylglycerophosphoethanolamine (Delta 9)[c]
% PEPALM_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoylglycerophosphoethanolamine[c]
% PE2LINL_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Linoleoylglycerophosphoethanolamine[c]
% PEAR_HSPLA2       H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Arachidonoyl-Sn-Glycero-3-Phosphoethanolamine[c]
% PE203_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatrienoylglycerophosphoethanolamine (Delta 11, 14, 17), Lpe (20:3)[c]
% PE226_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexenoylglyceroethanolamine (Delta 4, 7, 10, 13, 16, 19), Lpe (22:6)[c]
% PE224_HSPLA2      H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosatetraenoyglycerophosphoethanolamine (22:4, Delta 7, 10, 13, 16)[c]
% PEDH203_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Dihomo-Linolenoylglycerophosphoethanolamine (20:3, Delta 8, 11, 14)[c]
% PEDH12_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Didecanoylglycerophosphoethanolamine (C12:0 Pe)[c]
% PEDH14_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Myristoylglycerophosphoethanolamine (C14:0 Pe)[c]
% PEDH161_HSPLA2    H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Hexadecenoylglycerophosphoethanolamine (C16:1 Pe, Delta 9)[c]
% PEDH13_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Tridecanoylglycerophosphoethanolamine (C13:0 Pe)[c]
% PEDH15_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Pentadecanoylglycerophosphoethanolamine (C15:0 Pe)[c]
% PEDH17_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Heptadecanoylglycerophosphoethanolamine (C17:0 Pe)[c]
% PELINL_HSPLA2     H2O[c] + PE-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Linoleoylglycerophosphoethanolamine (Delta 9,12)[c]
%  
% PLA2_2            H2O[c] + PC-LD pool[c] => 2-lysolecithin pool[c] + H+[c] + R Total 2 Position[c]
% PLA2_2e           H2O[s] + PC-LD pool[s] => H+[s] + 2-lysolecithin pool[s] + R Total 2 Position[s]
% PCHOLMYR_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Myristoylglycerophosphocholine[c]
% PCHOLOLE_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Oleoylglycerophosphocholine (Delta 9)[c]
% PCHOLPALME_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoleoylglycerophosphocholine (Delta 9)[c]
% PCHOLPALM_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Palmitoylglycerophosphocholine[c]
% PCHOLSTE_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Stearoylglycerophosphocholine[c]
% PCHOL2LINL_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Linoleoylglycerophosphocholine[c]
% PCHOL2OLE_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Oleoylglycerophosphocholine[c]
% PCHOL2PALM_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Palmitoylglycerophosphocholine[c]
% PCHOL2STE_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 2-Stearoylglycerophosphocholine[c]
% PCHOLN15_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Pentadecanoylglycerophosphocholine, Sn1-Lpc (15:0)[c]
% PCHOLAR_HSPLA2    H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Arachidonoyl-Glycero-3-Phosphocholine[c]
% PCHOLN183_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 9, 12, 15)[c]
% PCHOLN1836_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Octadeca-Trienoylglycerophosphocholine, Sn1-Lpc (18:3, Delta 6, 9, 12)[c]
% PCHOLN19_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Nonadecanoylglycerophosphocholine, Sn1-Lpc (19:0)[c]
% PCHOLN201_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosenoylglycerophosphocholine (Delta 11) ,Sn1-Lpc (20:1)[c]
% PCHOLN204_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatetraenoylglycerophosphocholine (Delta 8, 11, 14, 17), Sn1-Lpc (20:4)[c]
% PCHOLN205_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosapentenoylglycerophosphocholine (Delta 5, 8, 11, 14, 17), Sn1-Lpc (20:5)[c]
% PCHOLN224_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosatetraenoylglycerophosphocholine (Delta 7, 10, 13, 16), Sn1-Lpc (22:4)[c]
% PCHOLN225_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosapentenoylglycerophosphocholine (Delta 7, 10, 13, 16, 19), Sn1-Lpc (22:5)-W3[c]
% PCHOLN2254_HSPLA2 H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosapentenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16), Sn1-Lpc (22:5)-W6[c]
% PCHOLN226_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexenoylglycerophosphocholine (Delta 4, 7, 10, 13, 16, 19), Sn1-Lpc (22:6)[c]
% PCHOLN203_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Dihomo-Linolenoylglycerophosphocholine (20:3, Delta 8, 11, 14), Lysopc A C20:3[c]
% PCHOLN24_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Lignocericylglycerophosphocholine (24:0), Lysopc A C24[c]
% PCHOLN261_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C26:1 (Delta 5)[c]
% PCHOLN281_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C28:1 (Delta 5)[c]
% PCHOLN28_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + Lysopc A C28:0[c]
% PCHOLDOC_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Docosahexaenoylglycerophosphocholine[c]
% PCHOLDEIC_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosadienoylglycerophosphocholine (Delta 11,14)[c]
% PCHOLDET_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Eicosatrienoylglycerophosphocholine (Delta 11, 14, 17)[c]
% PCHOLHEP_HSPLA2   H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Heptadecanoylglycerophosphocholine[c]
% PCHOLLINL_HSPLA2  H2O[c] + PC-LD pool[c] => H+[c] + R Total 2 Position[c] + 1-Linoleoylglycerophosphocholine (Delta 9,12)[c]
%  
% LPS2e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + 1-acylglycerol-3P-LD-TG1 pool[s] + R Total[s]
% MAGLINL_HSe   H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Linoleoylglycerol[s]
% MAGOLE_HSe    H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Oleoylglycerol[s]
% LPS5e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Palmitoylglycerol[s]
% LPS6e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Stearoylglycerol[s]
% LPS7e         H2O[s] + 1,2-diacylglycerol-LD-TAG pool[s] => H+[s] + R Total[s] + 1-Arachidonoyl Glycerol[s]

% load list of the above reactions to constrain
constrain_rxns = importdata('ComplementaryScripts/modelCuration/variable_mass_rxns_to_constrain.txt');
constrain_ind = ismember(ihuman.rxns,constrain_rxns);

% set upper and lower bounds to zero
ihuman.ub(constrain_ind) = 0;
ihuman.lb(constrain_ind) = 0;
ihuman.rev(constrain_ind) = 0;  % update reversibility





