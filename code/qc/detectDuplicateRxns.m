function repetitiveRxns=detectRepetitiveRxns(model,distReverse,metsToIgnore,distCoefficient)
% detectRepetitiveRxns
%   Detect repetitive reactions in a model
%
%   model                  a model structure
%   distReverse            distinguish reactions with same metabolites
%                          but different reversibility as different
%                          reactions (opt, default true)
%   metsToIgnore           either a cell array of metabolite IDs, a vector of 
%                          metabolite indexes to remove, or a logical vector with
%                          the same number of metabolites in the model
%   distCoefficient        disregard coefficient numbers in reactions by treating
%                          (opt, default false)
%
%
%   NOTE: This function detects duplicaitons only based on the stoicheometrix netwrok
%         regardless of other constraints
%
%   Usage: repetitiveRxns=detectRepetitiveRxns(model,distReverse,metsToIgnore,distCoefficient)
%


repetitiveRxns={};

if nargin<2
    distReverse=true;
end

if nargin<3
		skipMets=false;
else
	  if isempty(metsToIgnore)
	      skipMets=false;
	  else    
        indexesToIgnore=getIndexes(model,metsToIgnore,'mets');
        skipMets=true;
    end
end

if nargin<4
    distCoefficient=false;
end

%Construct equations for output
Eqns=constructEquations(model);

%If there are mets to ignore
if skipMets
    model.S(indexesToIgnore,:)=[];
end

%If disregard coefficient number
if distCoefficient
    model.S=spones(model.S);
end

%Transpose the matrix
if distReverse
    T=[model.S; model.rev']';
else
    T=model.S';
end

%Count occurrence of unique rxns
[~, I, J]=unique(T,'rows','first');
duplicateRxns=setdiff(1:numel(model.rxns),I);
uniqueRxns=I(J(duplicateRxns));

%Initialize cell array of repetitive rxn groups
repetitiveRxns.group=cell(numel(model.rxns),1);
repetitiveRxns.group(:)={''};
repetitiveRxns.equation=cell(numel(model.rxns),1);
repetitiveRxns.equation(:)={''};
repetitiveRxns.grRule=cell(numel(model.rxns),1);
repetitiveRxns.grRule(:)={''};

%Find repetitive rxn groups
for i=1:numel(duplicateRxns)
    %Generate cell array of repetitive reactions
    if ~isequal(duplicateRxns(i),uniqueRxns(i))
        %repetitiveRxns{uniqueRxns(i)}=[model.rxns{uniqueRxns(i)};model.rxns{duplicateRxns(i)];
        if isempty(repetitiveRxns.group{uniqueRxns(i)})
            %repetitiveRxns{uniqueRxns(i)}=[model.rxns{uniqueRxns(i)};model.rxns{duplicateRxns(i)}];
            repetitiveRxns.group{uniqueRxns(i)}{1,1}=model.rxns{uniqueRxns(i)};
            repetitiveRxns.group{uniqueRxns(i)}{2,1}=model.rxns{duplicateRxns(i)};
            repetitiveRxns.equation{uniqueRxns(i)}{1,1}=Eqns{uniqueRxns(i)};
            repetitiveRxns.equation{uniqueRxns(i)}{2,1}=Eqns{duplicateRxns(i)};
            repetitiveRxns.grRule{uniqueRxns(i)}{1,1}=model.grRules{uniqueRxns(i)};
            repetitiveRxns.grRule{uniqueRxns(i)}{2,1}=model.grRules{duplicateRxns(i)};
        else
            repetitiveRxns.group{uniqueRxns(i)}=[repetitiveRxns.group{uniqueRxns(i)};model.rxns{duplicateRxns(i)}];
            repetitiveRxns.equation{uniqueRxns(i)}=[repetitiveRxns.equation{uniqueRxns(i)};Eqns{duplicateRxns(i)}];
            repetitiveRxns.grRule{uniqueRxns(i)}=[repetitiveRxns.grRule{uniqueRxns(i)};model.grRules{duplicateRxns(i)}];
        end
    end
end

repetitiveRxns.group(strcmp('',repetitiveRxns.group))=[];
repetitiveRxns.equation(strcmp('',repetitiveRxns.equation))=[];
repetitiveRxns.grRule(strcmp('',repetitiveRxns.grRule))=[];

if isempty(repetitiveRxns.group)
    fprintf(['NO repetitive reactions found!\n']);
    %fprintf(['NO REPETITIVE REACTIONS FOUND\n']);
else
    fprintf([num2str(numel(repetitiveRxns.group)) ' groups of repetitive reactions found\n']);
end

end
