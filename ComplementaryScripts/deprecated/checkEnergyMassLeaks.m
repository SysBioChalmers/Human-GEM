function [outcome,results] = checkEnergyMassLeaks(model,printReport)
% checkEnergyMassLeaks
%   Verify that a model is properly balanced, such that mass and energy
%   cannot be created from nothing, or in unrealistic amounts. This is not
%   meant to be an exhaustive evaluation of a model's balance status, but a
%   quick check of some key features.
%
% Input:
%   model           model structure (must be HMR-derived, as the reaction
%                   names used are specific to HMR). 
%
%   printReport     true if a report should be printed to the screen
%                   (opt, default true)
%
% Output:
%   outcome         overall outcome; if all checks are passed, then outcome
%                   is 'PASS', otherwise it is 'FAIL'.
%
%   results         results structure containing information from each of
%                   the different tests
%       states      list of model states to examine
%                        'closed'  all exchange rxns closed
%                   'export_only'  exchange rxns allow only export
%                           'glc'  exchange rxns allow only export, and
%                                  glucose can be imported to max flux of 1
%                        'glc_o2'  exchange rxns allow only export, and
%                                  glucose can be imported to max flux of 1
%                                  and oxygen imported to max flux of 1000
%       objectives  list of model objectives to evaluate
%                    'ATP_hydrolysis'  maximize flux of ATP hydrolysis rxn
%                      'ATP_synthase'  maximize flux of ATP synthase rxn
%                    'CO2_production'  maximize export of CO2
%                     'O2_production'  maximize export of O2
%       outcome     cell array of test outcomes ('TRUE' or 'FALSE') for
%                   each combination of states and objectives, where rows
%                   correspond to the different states, and columns to the
%                   different objectives.
%       objVal      matrix of absolute objective value for each combination
%                   of states and objectives, where rows correspond to the
%                   different states, and columns to the different
%                   objectives.
%
% Usage: [outcome,results] = checkEnergyMassLeaks(model,printReport);
%
% Jonathan Robinson, 2018-11-12
%


% handle input arguments
if nargin < 2
    printReport = true;
end

% remove boundary metabolites (x-compartment)
model = simplifyModel(model);

% identify all exchange/sink/demand rxns (rxns with only 1 metabolite)
exch_inds = find(sum(model.S ~= 0) == 1);
model.S(:,exch_inds) = -abs(model.S(:,exch_inds));  % make all exch rxns the same direction (met --> nothing)

% close all exchange reactions for now, and reset objective vector
model = setParam(model, 'eq', exch_inds, 0);
model.c(:) = 0;

% define model states to check
states = {'closed'        % all exchange rxns closed
          'export_only'   % export is allowed for all exchange rxns (no import)
          'glc'           % export is allowed, and glucose can be imported (lb = -1)
          'glc_o2'};      % export is allowed, and glucose and oxygen can be imported (lb = -1 and -1000, respectivel)

% define model objectives to test for each state
objectives = {'ATP_hydrolysis'  % hydrolysis of ATP
              'ATP_synthase'    % synthesis of ATP
              'CO2_production'  % production (export) of CO2
              'O2_production'}; % production (export) of O2
          
% initialize results structure
results.states = states;
results.objectives = objectives;
results.outcome = repmat({''},numel(states),numel(objectives));
results.objVal = NaN(numel(states),numel(objectives));

% iterate through a the different model states
for i = 1:length(states)
    
    switch states{i}
        case 'closed'
            % all exchange rxns closed
            printProg(printReport,'\nModel state: All exchange reactions closed.\n');
            model_state = model;
            obj_max = [0;0;0;0];  % all objectives should be zero for this case
        case 'export_only'
            % exchange rxns can only export metabolites
            printProg(printReport,'\nModel state: Exchange rxns can only export metabolites.\n');
            model_state = setParam(model,'ub',exch_inds,1000);
            obj_max = [0;0;0;0];  % all objectives should be zero for this case
        case 'glc'
            % exchange rxns can only export, but glucose can be imported at
            % max flux of -1
            printProg(printReport,'\nModel state: All export open; max glucose import = 1.\n');
            model_state = setParam(model,'ub',exch_inds,1000);
            model_state = setParam(model_state,'lb','HMR_9034',-1);
            obj_max = [10;10;3;0];  % relatively relaxed limits on max ATP production, and 3 CO2 per glc is possible
        case 'glc_o2'
            % exchange rxns can only export, but glucose and oxygen can be
            % imported at max flux of -1 and -1000, respectively
            printProg(printReport,'\nModel state: All export open; max glucose import = 1, max O2 import = 1000.\n');
            model_state = setParam(model,'ub',exch_inds,1000);
            model_state = setParam(model_state,'lb',{'HMR_9034','HMR_9048'},[-1,-1000]);
            obj_max = [40;40;6;0];  % relatively relaxed limits on max ATP production, and 6 CO2 per glc is possible
    end
    
    % iterate through the different objectives to test
    for j = 1:length(objectives)
        
        switch objectives{j}
            case 'ATP_hydrolysis'
                % maximize flux through ATP hydrolysis rxn
                printProg(printReport,'\tObjective: maximize ATP hydrolysis');
                model_obj = setParam(model_state,'obj','HMR_3964',1);
            case 'ATP_synthase'
                % maximize flux through ATP synthase rxn
                printProg(printReport,'\tObjective: maximize ATP production');
                model_obj = setParam(model_state,'obj','HMR_6916',1);
            case 'CO2_production'
                % maximize export of CO2
                printProg(printReport,'\tObjective: maximize export of CO2');
                model_obj = setParam(model_state,'obj','HMR_9058',1);
                if strcmp(states{i},'closed')
                    % allow export of CO2 if all export is blocked
                    model_obj = setParam(model_obj,'ub','HMR_9058',1000);
                end
            case 'O2_production'
                % maximize export of O2
                printProg(printReport,'\tObjective: maximize export of O2');
                model_obj = setParam(model_state,'obj','HMR_9048',1);
                if strcmp(states{i},'closed')
                    % allow export of O2 if all export is blocked
                    model_obj = setParam(model_obj,'ub','HMR_9048',1000);
                end
        end
    
        % run FBA
        sol = solveLP(model_obj);
        printProg(printReport,'\tObj. value = %.3f',abs(sol.f));
        results.objVal(i,j) = abs(sol.f);
        
        % check that objective is below max allowed value
        if abs(sol.f) > (obj_max(j) + 0.001)  % allow for some error
            printProg(printReport,'\t(FAIL)\n');
        else
            printProg(printReport,'\t(PASS)\n');
            results.outcome{i,j} = 'PASS';
        end
    
    end
end

% determine overall outcome
if all(ismember(results.outcome,'PASS'))
    outcome = 'PASS';
else
    outcome = 'FAIL';
end


end


function printProg(printReport,varargin)
% function for printing results, to avoid having to use an if-statement
% every time we want to print something
    if printReport
        fprintf(varargin{:});
    end
end


