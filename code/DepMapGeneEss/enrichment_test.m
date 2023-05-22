%This file is for convenience copied from the Human1 paper.
function [penr,pdep] = enrichment_test(pop,sample,successes)
%enrichment_test caculates p-value of enrichment of successes in a sample.
%
% [PENR,PDEP] = enrichment_test(POP,SAMPLE,SUCCESSES) evaluates the
% significance of enrichment (and depletion) of SUCCESSES in a SAMPLE drawn
% from population POP using the hypergeometric test.
%
%
%--------------------------------- INPUTS ---------------------------------
%
% pop           Vector of genes comprising the population from
%               which samples are drawn.
%
% sample        Vector of genes sampled from the population.
%
% successes     Vector of genes in the population that are defined as
%               "successes". 
%               
%               The function will test if these "success" genes are 
%               significantly depleted or enriched in the sample, given
%               that they were drawn from the population.
%
% EXAMPLE: If a metabolite set enrichment analysis (MSEA) is being
%          performed, then SAMPLE is the list of metabolites of interest,
%          and SUCCESSES is a metabolite set (e.g., TCA cycle metabolites).
%          POP is the list of all metabolites from which SAMPLE and
%          SUCCESSES are drawn.
%
%--------------------------------- OUTPUTS --------------------------------
%
% penr          p-value associated with enrichment of successes in sample.
%
% pdep          p-value associated with depletion of successes in sample.
%
%
% ***** WARNING: P-VALUES ARE NOT ADJUSTED FOR FALSE DISCOVERY RATE! *****
%

x = numel(intersect(successes,sample));  % calc # of successes in sample
m = numel(pop);  % calc size of population
k = numel(intersect(successes,pop));
n = numel(sample);

penr = hygecdf(x-1,m,k,n,'upper');
pdep = hygecdf(x,m,k,n);
