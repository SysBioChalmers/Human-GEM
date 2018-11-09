%
% FILE NAME:    constrainVariableMassReactions.m
% 
% DATE CREATED: 2018-09-28
%     MODIFIED: 2018-11-06
% 
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: Script to constrain all reactions that involve an identical set
%          of metabolites except for one, and that one different metabolite
%          does not have the same mass in each reaction. For example, the
%          following reactions were identified as a set of "mass variable
%          reactions":
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


%% List of reactions to be constrained

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

% ARTFR13   myristoyl-CoA[c] => 0.875 R Group 1 Coenzyme A[c]
% ARTFR202  2 FADH2[m] + H+[m] + linolenoyl-CoA[c] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR203  2 FADH2[m] + gamma-linolenoyl-CoA[c] + H+[m] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR204  (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR205  dihomo-gamma-linolenoyl-CoA[c] + 2 FADH2[m] + H+[m] + NADPH[m] => 2 FAD[m] + NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR206  arachidonyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR207  eicosanoyl-CoA[c] => 1.25 R Group 2 Coenzyme A[c]
% ARTFR208  (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.25 R Group 2 Coenzyme A[c]
% ARTFR209  (7Z,10Z,13Z,16Z)-docosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR210  FADH2[m] + H+[m] + linoleoyl-CoA[c] + NADPH[m] => FAD[m] + NADP+[m] + 1.125 R Group 2 Coenzyme A[c]
% ARTFR211  (7Z,10Z,13Z,16Z,19Z)-docosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR212  (4Z,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR213  (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + 3 FADH2[m] + 3 H+[m] + 3 NADPH[m] => 3 FAD[m] + 3 NADP+[m] + 1.375 R Group 2 Coenzyme A[c]
% ARTFR31   stearoyl-CoA[c] => 1.125 R Group 3 Coenzyme A[c]
% ARTFR32   FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR33   FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR34   (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 3 Coenzyme A[c]
% ARTFR42   FADH2[m] + oleoyl-CoA[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR43   FADH2[m] + 11-Octadecenoyl Coenzyme A[c] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR44   (6Z,9Z)-octadecadienoyl-CoA[c] + 2 FADH2[m] => 2 FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR45   (15Z)-tetracosenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.5 R Group 4 Coenzyme A[c]
% ARTFR46   (2E)-octadecenoyl-CoA[c] + FADH2[m] => FAD[m] + 1.125 R Group 4 Coenzyme A[c]
% ARTFR51   tetracosanoyl-CoA[c] => 1.5 R Group 5 Coenzyme A[c]
% ARTFR52   hexacosanoyl-CoA[c] => 1.625 R Group 5 Coenzyme A[c]
% ARTFR53   (8Z,11Z,14Z,17Z)-eicosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.25 R Group 5 Coenzyme A[c]
% ARTFR54   (6Z,9Z,12Z,15Z,18Z)-tetracosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR55   (9Z,12Z,15Z,18Z,21Z)-tetracosapentaenoyl-CoA[c] + 3 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 3 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR56   (9Z,12Z,15Z,18Z)-tetracosatetraenoyl-CoA[c] + 2 FADH2[m] + 2 H+[m] + 2 NADPH[m] => 2 FAD[m] + 2 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
% ARTFR57   (6Z,9Z,12Z,15Z,18Z,21Z)-tetracosahexaenoyl-CoA[c] + 3 FADH2[m] + 3 H+[m] + 3 NADPH[m] => 3 FAD[m] + 3 NADP+[m] + 1.5 R Group 5 Coenzyme A[c]
%
% TAG_HSad      2 H2O[c] + 2 linoleoyl-CoA[c] + 2 oleoyl-CoA[c] + palmitoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] + stearoyl-CoA[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
% TAG_HSad_NE   2 H2O[c] + myristoyl-CoA[c] + 2 oleoyl-CoA[c] + palmitoleoyl-CoA[c] + palmitoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] + stearoyl-CoA[c] => 6 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
% TAG_HSad_E    (4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-CoA[c] + (5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c] + arachidonyl-CoA[c] + 2 H2O[c] + linolenoyl-CoA[c] + linoleoyl-CoA[c] + 2 sn-glycerol-3-phosphate[c] => 5 CoA[c] + 2 Pi[c] + 2 TAG-VLDL pool[c]
%
% This reaction was found to exhibit variable mass with other reactions
% involving cholesterol and cholesterol-ester pools in the model:
% CHOLESTle     cholesterol-ester pool[s] + H2O[s] => cholesterol[s] + H+[s] + R Total[s]
%

%% Constrain reactions

% load model if it does not yet exist
if ~exist('ihuman','var')
    load('humanGEM.mat');  % version 0.5.2
end
ihuman_orig = ihuman;

% load list of the above reactions to constrain (stored in txt file)
constrain_rxns = importdata('../../ComplementaryData/modelCuration/variable_mass_rxns_to_constrain.tsv');
constrain_ind = ismember(ihuman.rxns,constrain_rxns);
if any(~ismember(constrain_rxns,ihuman.rxns))
    error('Some reactions were not found in the model.');
end

% generate rxnNotes array
rxnNotes = [constrain_rxns, repmat({'reaction treats the mass of one or more of its metabolites in an inconsistent manner, resulting in mass imbalances; rxn should therefore be constrained until imbalances can be addressed, otherwise DELETED'},length(constrain_rxns),1)];

% set upper and lower bounds to zero
ihuman.ub(constrain_ind) = 0;
ihuman.lb(constrain_ind) = 0;
ihuman.rev(constrain_ind) = 0;  % update reversibility

% document reaction changes
rxnChanges = docRxnChanges(ihuman_orig,ihuman,rxnNotes);
writeRxnChanges(rxnChanges,'constrainVariableMassReactions_rxnChanges',true);

% remove intermediate variables
clearvars -except ihuman

% save new model file
save('../../ModelFiles/mat/humanGEM.mat','ihuman');
movefile('constrainVariableMassReactions_rxnChanges.tsv','../../ComplementaryData/modelCuration/');

