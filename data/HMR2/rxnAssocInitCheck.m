%
%   FILE NAME:    rxnAssocInitCheck.m
% 
%   PURPOSE: Initial consistence check for associated exteranl reaction
%            identifiers in HMR2. This task was triggered by the discovery of
%            mistakenly-associated HepatoNet1 reactions to BiGG. This script
%            focus on curating the identifiers derived from EHMN and HepatoNet1
%            using the orgianl publications/files and database investigation.
%


% 1. Load HMR model ver 2.0.1
load('HMRdatabase2_01.mat');


% 2. Check consistence of EHMN reaction association in HMR2
% This check is based on that EHMN and KEGG rxns are closely associated 

for i=1:numel(ihuman.rxns)
	% Detect if their identifiers are consistent
	if ~isempty(ihuman.rxnEHMNID{i})

		% A. Check EHMN identifiers without KEGG association (14)
		if isempty(ihuman.rxnKEGGID{i})
			if (regexp(ihuman.rxnEHMNID{i},'R\d+'))
				% Find out the KEGG
				%disp([num2str(i) ':[' ihuman.rxnEHMNID{i} '-' ihuman.rxnKEGGID{i} ']']);
					%556:[R00707M-]
					%735:[R03102M-]
					%1677:[R07057-]
					%1712:[R07770-]
					%1713:[R07771-]
					%1714:[R07768-]
					%1715:[R07769-]
					%1716:[R07766-]
					%1717:[R07767-]
					%1718:[R08550-]
					%2221:[R00521X-]
					%2222:[R00521C-]
					%2925:[R02208-]
					%2926:[R08957-]
					%Manually fix them later
			else
				%Ignore the (1186) cases that have prefix RE (1092), RT (91) and RN (3) 
				%disp(['[' ihuman.rxnEHMNID{i} '-' ihuman.rxnKEGGID{i} ']']);
			end

		% Check EHMN identifiers with KEGG association
		else
			% B. If they are identical, then remove EHMN ones that appeared strange (12)
			if isequal(ihuman.rxnKEGGID{i},ihuman.rxnEHMNID{i})
				disp(['[' ihuman.rxnEHMNID{i} '-' ihuman.rxnKEGGID{i} ']']);
				%[R01252-R01252]
				%[R01252-R01252]
				%[R04545-R04545]
				%[R08379-R08379]
				%[R02124-R02124]
				%[R08379-R08379]
				%[R02124-R02124]
				%[R08387-R08387]
				%[R08388-R08388]
				%[R08389-R08389]
				%[R02123-R02123]
				%[R02123-R02123]
				ihuman.rxnEHMNID{i}='';
			% If they are identical after removing last comp character of EHMN id
			elseif isequal(ihuman.rxnKEGGID{i},regexprep(ihuman.rxnEHMNID{i},'\w$',''))
				%These are normal cases and ignore them
				%disp(['[' ihuman.rxnEHMNID{i} '-' ihuman.rxnKEGGID{i} ']']);
			% C. There are (51) probalmatic cases and require manual curation
			else
				%disp([num2str(i) ':[' ihuman.rxnEHMNID{i} '-' ihuman.rxnKEGGID{i} '-' ihuman.rxnHepatoNET1ID{i} ']']);
			end
		end
		
	end
end

% A. Manually fix EHMN identifiers without KEGG association (14)
%leave it                         %556:[R00707M-]
ihuman.rxnKEGGID{735}='R03102';   %735:[R03102M-]
ihuman.rxnKEGGID{2221}='R00521';  %2221:[R00521X-]
ihuman.rxnKEGGID{2222}='R00521';  %2222:[R00521C-]
ihuman.rxnKEGGID{1677}='R07057';  %1677:[R07057-]
ihuman.rxnEHMNID{1677}='R07057C'; %1677:[R07057-]
ihuman.rxnKEGGID{1712}='R07770';  %1712:[R07770-]
ihuman.rxnEHMNID{1712}='R07770C'; %1712:[R07770-]
ihuman.rxnKEGGID{1713}='R07771';  %1713:[R07771-]
ihuman.rxnEHMNID{1713}='R07771C'; %1713:[R07771-]
ihuman.rxnKEGGID{1714}='R07768';  %1714:[R07768-]
ihuman.rxnEHMNID{1714}='R07768C'; %1714:[R07768-]
ihuman.rxnKEGGID{1715}='R07769';  %1715:[R07769-]
ihuman.rxnEHMNID{1715}='R07769C'; %1715:[R07769-]
ihuman.rxnKEGGID{1716}='R07766';  %1716:[R07766-]
ihuman.rxnEHMNID{1716}='R07766C'; %1716:[R07766-]
ihuman.rxnKEGGID{1717}='R07767';  %1717:[R07767-]
ihuman.rxnEHMNID{1717}='R07767C'; %1717:[R07767-]
ihuman.rxnKEGGID{1718}='R08550';  %1718:[R08550-]
ihuman.rxnEHMNID{1718}='R08550C'; %1718:[R08550-]
ihuman.rxnKEGGID{2925}='R02208';  %2925:[R02208-]
ihuman.rxnEHMNID{2925}='R02208C'; %2925:[R02208-]
ihuman.rxnKEGGID{2926}='R08957';  %2926:[R08957-]
ihuman.rxnEHMNID{2926}='R08957C'; %2926:[R08957-]

% C. Manually fix the 51 probalmatic cases with conflicting EHMN and KEGG ids
ihuman.rxnEHMNID{42}='';          %42:[RE0446-R00028-r0014]
ihuman.rxnKEGGID{44}='R00010';    %44:[R00010C-R06103-r0785]
ihuman.rxnEHMNID{165}='';         %165:[RE2649C-R00925-r0218
ihuman.rxnBiGGID{165}='ACS2';     %165:[RE2649C-R00925]
ihuman.rxnKEGGID{305}='';         %305:[RE2813C-R02425]
ihuman.rxnEHMNID{729}='R01939C';  %729:[R01904C-R01939-r0450]
%leave it                         %809:[RE3326M-R02662-r0561]
%leave it                         %892:[RE1927C-R04084-]  Have different MNXref assoc
%leave it                         %966:[RE1915C-R08848-]      
%leave it                         %1010:[RE0159C-R03106-]     
%leave it                         %1679:[RE3409C-R05718-]     
%leave it                         %1928:[RE3038C-R07032-]  use KEGG assoc:R07032
%leave it                         %1929:[RE3038N-R07032-]  use KEGG assoc:R07032
%leave it                         %1930:[RE3038R-R07032-]  use KEGG assoc:R07032
%leave it                         %1931:[RE3038X-R07032-]  use KEGG assoc:R07032
%leave it                         %1932:[RE3040C-R07031-]  use KEGG assoc:R07031
%leave it                         %1933:[RE3040N-R07031-]  use KEGG assoc:R07031
%leave it                         %1934:[RE3040R-R07031-]  use KEGG assoc:R07031
%leave it                         %1935:[RE3040X-R07031-]  use KEGG assoc:R07031
%leave it                         %1945:[RE3520C-R07039-]  use KEGG assoc:R07039
%leave it                         %1987:[RE3010C-R03863-]  use KEGG assoc:R03863
%leave it                         %1988:[RE3010M-R03863-]  use KEGG assoc:R03863
%leave it                         %1989:[RE3470M-R03864-]  use KEGG assoc:R03864
%leave it                         %1990:[RE3010R-R03863-]  use KEGG assoc:R03863
%leave it                         %1991:[RE3470X-R03864-]  use KEGG assoc:R03864
%leave it                         %1992:[RE3010X-R03863-]  use KEGG assoc:R03863
%leave it                         %2125:[RE3550X-R04256-]  use KEGG assoc:R04256
%leave it                         %2192:[RE0579C-R08185-]  use KEGG assoc:R08185
%leave it                         %2643:[RE3075C-R03631-]  use KEGG assoc:R03631
%leave it                         %2645:[RE3075X-R03631-]  use KEGG assoc:R03631
ihuman.rxnEHMNID{2798}='';        %2798:[RE0511M-R04100-r0657]  use KEGG assoc:R04100
%leave it                         %2901:[RE1100C-R08941-]  use KEGG assoc:R08941
ihuman.rxnKEGGID{3032}='';        %3032:[RE2410C-R03724-r0632]  KEGG assoc was wrong
%leave it                         %3053:[RE3136C-R03353-r1381]  use KEGG assoc:R03353
ihuman.rxnEHMNID{3486}='';        %3486:[R00848C-R00849-r0205]  complicated case, inaccurate KEGG assoc and wrong EHMN assoc to HepatoNet1
%leave it                         %3537:[RE0066C-R03424-]  the two are consistent
%leave it                         %3538:[RE3511C-R01320-]  the two are consistent
ihuman.rxnEHMNID{3671}='';        %3671:[RE2679-R04018-]
%leave it                         %3735:[RE2078R-R02583-]  use KEGG assoc:R02583
%leave it                         %3740:[RE2799C-R02264-]  use KEGG assoc:R02264, balance adjusted in HMR2
%leave it                         %3745:[RE2079R-R02801-]  use KEGG assoc:R02801
%leave it                         %3746:[RE2079R-R02801-]  use KEGG assoc:R02801
%leave it                         %3793:[RE3556C-R04565-]  the two are consistent
ihuman.rxnKEGGID{4105}='';        %4105:[RE3247X-R04592-r0706]  KEGG assoc was wrong
ihuman.rxnKEGGID{4113}='R07296';  %4112:[RE1834C-R07296-r0794]  use KEGG assoc:R07296
ihuman.rxnKEGGID{4114}='R07296';  %4112:[RE1834C-R07296-r0794]  use KEGG assoc:R07296
%leave it                         %4433:[RE2974C-R05802-]  use KEGG assoc:R05802
%leave it                         %4435:[RE3272N-R03361-]  4435 and 4436 are two identical rxns
%leave it                         %4436:[RE3272N-R03361-]  with opposite direction, strange
ihuman.rxnEHMNID{4443}='RE2972M'; %4443:[RE1447M-R05803-]  correct association error
%leave it                         %4751:[RE2426C-R03629-]  Conflicting! KEGG rxn has more consistent mets but different equation
ihuman.rxnEHMNID{4753}='';        %4753:[RE2439C-R03628-]  use KEGG assoc:R03628
ihuman.rxnKEGGID{5022}='R03538';  %5022:[RE1860C-R03422-]  use EHMN assoc:RE1860C and replace KEGG assoc with R03538

% Many above problems appeared as typos
save('HMRdatabase2_02.mat','ihuman');  %===2018-02-03


% 3. Check consistence of HepatoNet1 reaction association in HMR2

% Load HepatoNET1 reaction info for investigation
% The reactions data was obtained from the supplementary information
% of the publication that was downloaded from the journal website
T=readtable('inline-supplementary-material-4.txt','Delimiter','tab');
HepatoNet1=table2struct(T,'ToScalar',true);
% Remove a dash character found as prefix in some Recon1 identifiers
HepatoNet1.Recon1=regexprep(HepatoNet1.Recon1,'^\-','');
save('HepatoNet1.mat','HepatoNet1');  %2018-01-25 

% Load Recon1 reaction identifiers for investigation
[~, textData]=xlsread('Recon1_rxns.xlsx','Sheet1');
Recon1.rxns=textData(3:end,1);
save('Recon1_rxns.mat','Recon1_rxns');  %2018-01-25 

% Loop through all rxns
%count=0;
for i=1:numel(ihuman.rxns)
	% Go through existing HepatoNet1 reaction assocations
	if ~isempty(ihuman.rxnHepatoNET1ID{i})
		[a, b]=ismember(ihuman.rxnHepatoNET1ID{i},HepatoNet1.r_ID);
		if a
			% Check if KEGG associations are consistent between HMR2 and HepatoNet1
			if ~isempty(HepatoNet1.KEGG{b})  %HepatoNet1 has KEGG assoc
				% A. These HMR2 rxns have no KEGG assoc (21)
				if isempty(ihuman.rxnKEGGID{i})
					%disp([num2str(i) ': ' HepatoNet1.KEGG{b}]);
					%207: R01847   !This rxn has been removed by KEGG, what to do?!
					%556: R00707   The KEGG assocs of HepatoNet1 was wrong, also EHMN assoc
					%3030: R04804  Add this KEGG assoc
					%3032: R03724  leave out this KEGG assoc
					%4080: R04817  Add this KEGG assoc
					%4092: R04507  Add this KEGG assoc, the BiGG assoc was wrong!
					%4098: R04826  The HepatoNet1 assoc was wrong, remove it for 4098, 4099, 4100
					%4099: R04826  also found wrong HepatoNet1 assoc for 4101 and 4102
					%4100: R04826  Fix assoc to 4101 and 4102 through using REACTOME assoc
					%4105: R04592  wrong KEGG assoc, already removed based on EHMN; find same rxn for 4103 then fix it
					%4115: R04580  wrong HepatoNet1 assoc, remove it
					%4116: R04580  wrong HepatoNet1 assoc, remove it
					%4124: R04506  wrong HepatoNet1 assoc, remove it
					%4125: R04506  wrong HepatoNet1 assoc, remove it
					%4127: R04507  wrong HepatoNet1 assoc, remove it
					%4133: R04546  wrong HepatoNet1 assoc, remove it
					%4141: R04823  wrong HepatoNet1 assoc, remove it
					%4142: R04825  wrong HepatoNet1 assoc, remove it
					%4145: R03506  wrong HepatoNet1 assoc, remove it
					%4410: R01623  Add this KEGG assoc
					%6777: R04806  wrong HepatoNet1 assoc, remove it
					%Use HepatoNet1 KEGG assoc to fill these, with manual curation
				% B. These HMR2 rxns have conflicting KEGG assoc with HepatoNet1 (5)
				elseif ~isequal(ihuman.rxnKEGGID{i},HepatoNet1.KEGG{b})  % conflicting cases
					%disp([num2str(i) ': [' ihuman.rxnKEGGID{i} '-' HepatoNet1.KEGG{b} ']']);
					%802: [R01214-R01090]   leave out this wrong KEGG assoc to HepatoNet1
					%3486: [R00849-R00848]  leave out this wrong KEGG assoc to HepatoNet1, EHMN assoc already removed
					%4084: [R04807-R04805]  wrong HepatoNet1 assoc, replace with r0742
					%4106: [R04813-R04807]  wrong HepatoNet1 assoc, replace with r0744
					%4120: [R04817-R04818]  wrong HepatoNet1 assoc, replace with r0745
					%Manually check above cases
				end
			end

			% Check if Recon1 associations are consistent between HMR2 and HepatoNet1 (10)
			if ~isempty(HepatoNet1.Recon1{b})  %HepatoNet1 has BiGG assoc
				%C. HMR2 rxns have no Recon1 associaiton
				if isempty(ihuman.rxnBiGGID{i})
					%disp([num2str(i) ': ' HepatoNet1.Recon1{b}]);
					%4080: AKR1D    Add this Recon1 assoc
					%4098: P4508B11r   the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4099: P4508B11r   the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4100: P4508B11r   the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4101: XOLDIOLONEt the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4115: VLCSr       the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4116: VLCSr       the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4127: VLCS2r      the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4141: AKR1D2      the HepatoNet1 assoc was wrong, so ignore its Recon1 assoc
					%4143: XOLTRIOLtm  wrong HepatoNet1 assoc, remove it
					%Use HepatoNet1 Recon1 assoc to fill these, after manual curation
					
				elseif isequal(ihuman.rxnBiGGID{i},HepatoNet1.Recon1{b})
				% If Recon1 assoc in HepatoNet1 is consistent with HMR2 (600)
				% These are normal cases and ignore them
					%count=count+1;  % 600 cases
				
				% When Recon1 assoc in HepatoNet1 is not consistent with HMR2
				else
					%disp([num2str(i) ': [' ihuman.rxnBiGGID{i} '-' HepatoNet1.Recon1{b} ']']);
					%If the Recon1 assoc of these HepatoNet1 rxns are found in Recon1
					if ~ismember(HepatoNet1.Recon1{b},Recon1_rxns)
						%disp([num2str(i) ': [' ihuman.rxnBiGGID{i} '-' HepatoNet1.Recon1{b} ']']);
						%No hits, this means Recon1 assoc of these HepatoNet1 rxns are all consistent
						
					else
					%Then check Recon1 assoc of these HMR rnxs here
						if ismember(ihuman.rxnBiGGID{i},Recon1_rxns)
							%D. Some are also found in Recon1 (8)
							%disp([num2str(i) ': [' ihuman.rxnBiGGID{i} '-' HepatoNet1.Recon1{b} ']']);
							%184: [ACACT1r-ACACT1]         updated BiGG id, leave it
							%4120: [AKR1D-AKR1C41]         updated BiGG id, leave it
							%4581: [MTHFD2-MTHFD2m]        replace with HepatoNet1 Recon1 assoc
							%6777: [XOLTRIOLtm-P45027A14m] updated BiGG id, leave it
							%7595: [H2Oter-H2Otg]          updated BiGG id, leave it
							%7598: [CO2ter-COAtr]          replace with HepatoNet1 Recon1 assoc
							%7605: [GLCter-GLCtg]          updated BiGG id, leave it
							%7613: [PIter-PItg]            updated BiGG id, leave it
							%Manually check above cases
						else
							%E. The BiGG/Recon1 assocs of these HMR rxns should be wrong and listed in the end (119)
							disp([num2str(i) ': [' ihuman.rxnBiGGID{i} '-' HepatoNet1.Recon1{b} ']']);
							%they are corrected by the Recon1 assoc of corresponding HepatoNet1 and listed below
							ihuman.rxnBiGGID{i}=HepatoNet1.Recon1{b};
						end
					end
				end
			end
			
		% F. Detect identifiers that were not found in HepatoNet1 model and were
		% found are KEGG rxn ids. They (R00736, R02384, R02382, R02695, R02697)
		% were fixed by moving to KEGG ids after manual curation (5)
		else
			disp([num2str(i) ': ' ihuman.rxnHepatoNET1ID{i}]);
			%913: R00736
      %914: R02384
      %915: R02382
      %916: R02695
      %917: R02697
			ihuman.rxnKEGGID{i}=ihuman.rxnHepatoNET1ID{i};
			ihuman.rxnHepatoNET1ID{i}='';
		end
	end
end

%Manual curations:
% A. fix HMR2 rxns have no KEGG assoc, fix them manually (21)
%Remove this reaction?!                  %207: R01847   !This rxn has been removed by KEGG, what to do?!
ihuman.rxnKEGGID{556}='R03314';          %556: R00707   The KEGG assocs of HepatoNet1 was wrong
ihuman.rxnEHMNID{556}='';                %556: R00707   The EHMN assoc was also wrong
ihuman.rxnKEGGID{3030}='R04804';         %3030: R04804  Add this KEGG assoc
%leave it                                %3032: R03724  leave out this KEGG assoc
ihuman.rxnKEGGID{4080}='R04817';         %4080: R04817  Add this KEGG assoc
ihuman.rxnKEGGID{4092}='R04507';         %4092: R04507  Add this KEGG assoc, the BiGG assoc was wrong!
ihuman.rxnHepatoNET1ID{4098}='';         %4098: R04826  The HepatoNet1 assoc for 4098 was wrong, remove it
ihuman.rxnHepatoNET1ID{4099}='';         %4099: R04826  The HepatoNet1 assoc for 4099 was wrong, remove it
ihuman.rxnHepatoNET1ID{4100}='';         %4100: R04826  The HepatoNet1 assoc for 4100 was wrong, remove it
ihuman.rxnHepatoNET1ID{4101}='';         %4099: R04826  wrong HepatoNet1 assoc for 4101 and 4102
ihuman.rxnREACTOMEID{4102}='REACT_9998'; %4100: R04826  Fix assoc to 4101 and 4102 through REACTOME assoc
ihuman.rxnHepatoNET1ID{4103}='r0706';    %4103: R04592  Fix assoc to 4103 through REACTOME and HepatoNet1
ihuman.rxnREACTOMEID{4103}='REACT_10074';%4103: R04592  Fix assoc to 4103 through REACTOME and HepatoNet1
ihuman.rxnHepatoNET1ID{4115}='';         %4115: R04580  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4116}='';         %4116: R04580  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4124}='';         %4124: R04506  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4125}='';         %4125: R04506  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4127}='';         %4127: R04507  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4133}='';         %4133: R04546  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4141}='';         %4141: R04823  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4142}='';         %4142: R04825  wrong HepatoNet1 assoc, remove it
ihuman.rxnHepatoNET1ID{4145}='';         %4145: R03506  wrong HepatoNet1 assoc, remove it
ihuman.rxnKEGGID{4410}='R01623';         %4410: R01623  Add this KEGG assoc
ihuman.rxnHepatoNET1ID{6777}='';         %6777: R04806  wrong HepatoNet1 assoc, remove it

% B. These HMR2 rxns have conflicting KEGG assoc with HepatoNet1, fix them manually (5)
ihuman.rxnHepatoNET1ID{4084}='r0742';    %4084: [R04807-R04805]  wrong HepatoNet1 assoc, replace with r0742
ihuman.rxnHepatoNET1ID{4106}='r0744';    %4106: [R04813-R04807]  wrong HepatoNet1 assoc, replace with r0744
ihuman.rxnHepatoNET1ID{4120}='r0745';    %4120: [R04817-R04818]  wrong HepatoNet1 assoc, replace with r0745
ihuman.rxnKEGGID{4119}='R03718';         %4119: R03720  The KEGG assocs of HepatoNet1 was wrong and replace with R03718
ihuman.rxnEHMNID{4119}='';               %4119: R03720  The EHMN assocs of HepatoNet1 was also wrong

%C. HMR2 rxns have no Recon1 associaiton, fix them manually (10)
ihuman.rxnBiGGID{4080}='AKR1D';          %4080: AKR1D  Add this Recon1 assoc
ihuman.rxnHepatoNET1ID{4143}='';         %4143: XOLTRIOLtm  wrong HepatoNet1 assoc, remove it

%D. Some are also found in Recon1 (8)
ihuman.rxnBiGGID{4581}='MTHFD2m';          %4581: [MTHFD2-MTHFD2m]  replace with HepatoNet1 Recon1 assoc
ihuman.rxnBiGGID{7598}='COAtr';            %7598: [CO2ter-COAtr]  replace with HepatoNet1 Recon1 assoc

%E. These HMR rxns BiGG/Recon1 assoc should be wrong and listed below (119)
%disp([num2str(i) ': [' ihuman.rxnBiGGID{i} '-' HepatoNet1.Recon1{b} ']']);
%they are corrected by the Recon1 assoc of corresponding HepatoNet1 and listed below
%ihuman.rxnBiGGID{i}=HepatoNet1.Recon1{b};
31: [G6PASEer-G6PPer]
139: [ME2r-ME2]
140: [ME2rm-ME2m]
143: [PPCKG-PEPCK]
144: [PPCKGm-PEPCKm]
167: [PPSm-ACCOALm]
192: [G6PDH2er-G6PDH2rer]
193: [G6PDHy-G6PDH2r]
201: [PGDH-GND]
205: [TAL-TALA]
241: [ADSL1r-ADSL1]
290: [PRAGS-PRAGSr]
295: [ADSL2r-ADSL2]
329: [TRDRr-TRDR]
365: [DHORD3m-ubq10-DHORD9]
369: [TMDSr-TMDS]
406: [DCMPDA2ir-DCMPDA]
413: [GARFTi-GARFT]
518: [ASPTA1m-ASPTAm]
519: [ASPTA1-ASPTA]
528: [GLUNc-GLUNm]
556: [G5SADsm-G5SADrm]
563: [G5SDm-G5SDym]
660: [PGCDr-PGCD]
661: [PSERTr-PSERT]
663: [GHMT-GHMT2r]
670: [GMTR-GNMT]
796: [ACACT10rm-ACACT10m]
797: [PPCOACrm-PPCOACm]
800: [MMMrm-MMMm]
904: [HPPDO1-34HPPOR]
905: [HGENDO-HGNTOR]
906: [MLACI-MACACI]
907: [FUMACA-FUMAC]
987: [ASPTA4-CYSTA]
1074: [MOBD2m-OIVD2m]
1080: [MOBD1m-OIVD1m]
1082: [MCCCm-MCCCrm]
1083: [MGCHm-MGCHrm]
1084: [MOBD3m-OIVD3m]
1086: [ECOAH3m-ECOAH9m]
1090: [HACD8m-HACD9m]
1357: [SUCOASGm-SUCOAS1m]
1358: [CITL2-ACITL]
1359: [SUCOASAm-SUCOASm]
1400: [SOD-SPODM]
1424: [FACOAL160-FACOAL160i]
1433: [FACOAL181n9-FACOAL181i]
1465: [FACOAL182n6-FACOAL1821]
1468: [FACOAL204n6-FACOAL204]
3008: [MEVK1p-MEVK1x]
3009: [PMEVKrp-PMEVKx]
3010: [IPDDIp-IPDDIx]
3013: [SQLEer-SQLEr]
3014: [LNSTLSer-LNSTLSr]
3030: [CHSTNIer-EBP1r]
3031: [LSTO1er-LSTO1r]
3032: [DHCR71er-DHCR71r]
3035: [LATHSTOxer-LSTO2r]
3036: [LATHSTOyer-LSTO2r]
3037: [DHCR72er-DHCR72r]
4073: [P4507A1er-CH25H]
4076: [P4508B11er-P4508B11r]
4079: [XOLDIOLONEter-XOLDIOLONEt]
4092: [THCLSTCtm-VLCS2p]
4238: [CDO-CYSO]
4373: [HANTHDOr-3HAO]
4391: [THD1im-THD1m]
4519: [DHPRx-DHPR]
4520: [FORTHFC-FTCD]
4525: [FMETDH-FTHFDH]
4526: [FTHFLr-FTHFL]
5212: [HDCAtb-HDCAtr]
5230: [OCDCEA9tb-OCDCEAtr]
5270: [OCDCTRA3tb-LNLNCAt]
5294: [OCDDEA6tb-LNLCt]
5300: [ECSTTEA6tb-ARACHDt2]
5317: [GLYCt5b-GLYCt]
5370: [O2tb-O2t]
5374: [CO2tb-CO2t]
5380: [ACt6bl-ACt2r]
5388: [CHOLtb-CHOLtu]
5407: [ETOHtb-ETOHt]
5418: [D3AIBt-D-3AIBt]
5420: [NH4tb-NH4t3r]
5421: [GLCt1b-GLCt2r]
5431: [ARGtrb-ARGtiDF]
5433: [LYStb-LYStiDF]
5435: [TYRtb-TYRt]
5440: [TRPtb-TRPt]
5441: [PHEtb-PHEtec]
5443: [LEUtb-LEUtec]
5446: [VALtb-VALtec]
5448: [ILEt5b-ILEtec]
5470: [ALAt4rb-ALAt4]
5472: [GLNt4rb-GLNt4]
5474: [ASNt4rb-ASNt4]
5482: [GLYt4rb-GLYt4]
5483: [PROt4rb-PROt4]
5484: [METt4rb-METt4]
5485: [THRt4rb-THRt4]
5489: [HISt4b-HISt4]
5493: [NAt7b-NAt3_1]
6044: [FE2tb-FE2t]
6053: [L-LACt2b-L-LACt2r]
6354: [SERt4rb-SERt4]
6355: [CYSt4rb-CYSt4]
6798: [GACm-ASPGLUm]
6825: [MALAKGtm-AKGMALtm]
6839: [O2trm-O2tm]
6843: [PYRtim-PYRt2m]
6850: [CITMALtm-CITtam]
6866: [Pitm-PIt2m]
6870: [GLNtrm-GLNtm]
6896: [LLACtm-L-LACtm]
6897: [FE2trm-FE2tm]
6940: [ATP/ADPtm-ATPtm]
7130: [PIt2p-PItx]
7601: [FORter-FORtr]


save('HMRdatabase2_02.mat','ihuman');  %===2018-02-07

%In sum, a total of 2031 reactions were investigated in this script.
%Inconsistent cases were detected in 1431 reactions, among which
%299 reactions were manually checked and corrected for their external
%identifiers while problems found in the other 1132 reactions were
%fixed automatically by this script.
