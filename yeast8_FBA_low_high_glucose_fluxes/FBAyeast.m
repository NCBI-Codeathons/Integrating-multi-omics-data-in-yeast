%FBA on yeast 

%software installation as described in https://opencobra.github.io/cobratoolbox/stable/installation.html

%FBA code modified from https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialFBA.html

%yeast model downloaded from https://github.com/SysBioChalmers/yeast-GEM

%Load the model
%modelyeast=loadYeastModel

%default contraints in the model:
%printConstraints(modelyeast,-1000,1000) 
%r_1714 lower bound = -1 : D-glucose exchange D-glucose[e] <=>  
%r_4046 lower and upper bound = 0.7 : 	non-growth associated maintenance reaction	ATP[c] + H2O[c] => ADP[c] + H+[c] + phosphate[c]


%by default the objective function is the biomass production reaction
%modelyeast.rxns(find(modelyeast.c)) %this is r_2111 which is biomass prod
%to set a different objective function reaction use:
%modelyeast = changeObjective(modelyeast,'r_xxxx');

%To perform FBA with maximization of the biomass reaction as the objective, enter:

FBAsolution_normal = optimizeCbModel(modelyeast,'max')

%to get the value of the objective function (biomass production)
%modelyeast.c'*FBAsolution_normal.v
%or equivalently
%FBAsolution_normal.f

%for the default model the value was 0.088

%now, reducing the glucose exchange:
lowglucoseuptake1 = -33.2*60/1000 % -1.992 mmol/gDW/hr millimoles per gram dry cell weight per hour
modelyeast_low = changeRxnBounds(modelyeast,'r_1714',lowglucoseuptake1,'l');
printConstraints(modelyeast_low,-1000,1000) 
FBAsolution_low = optimizeCbModel(modelyeast_low,'max')
%f was 0.072
fluxprofilelow=FBAsolution_low.v;
save('flux_profile_yeast8_lowglucose.txt','fluxprofilelow','-ascii')

%now, with highglucose exchange:
highglucoseuptake1 = -71.2*60/1000% -4.2720 mmol/gDW/hr
modelyeast_high = changeRxnBounds(modelyeast,'r_1714',highglucoseuptake1,'l');
printConstraints(modelyeast_high,-1000,1000) 
FBAsolution_high = optimizeCbModel(modelyeast_high,'max')
%f was 0.3859
fluxprofilehigh=FBAsolution_high.v;
save('flux_profile_yeast8_highglucose.txt','fluxprofilehigh','-ascii')

