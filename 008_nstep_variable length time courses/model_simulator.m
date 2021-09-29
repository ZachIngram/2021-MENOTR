function [sim_data_compiled]=model_simulator(params,indvar,flag)
% ALL SECTIONS THAT NEED TO BE MODIFIED ARE MARKED WITH ***
 
% Sections that require modification are:
    % Section 2
    % Section 3 
 
%% Section 1: Check Inputs  
% In this section of the code, the input parameter "params"  
% is checked to make sure they're in the right orientation before moving 
% forward to the next section.
[r,c]=size(params);
if r>c
    params=params';
end
 
 
%% Section 2: Define Independent Variable and Parameters*** 

%This section defines the parameters that will be used in your
    %model. Each parameter will be set equal to an element in a matrix
    %called params. Each parameter will be set equal to an element in a
    % matrix called params (Example: A = params(#); where A is a fitting parameter
    % and # is the position it occupies in the params matrix.) List
    %all of your parameters following the example given below. These will
    %be used in the model below. If you have local parameters in your fit
    %you will need a unique parameter for each data set for that parameter.
    %See A in example 1 if we had 3 data sets to fit locally for A.
% Example
%     k1 = params(1);
%     k2 = params(2);
%     m = params(3);
%     A_1 = params(4);
%     A_2 = params(4);
%     A_3 = params(4);
 
 
% Define your independent variable below.
% time = indvar;
for i = 1: length(flag)
    if flag(i) == 1
        fp_t1 = i;
    elseif flag(i) == 2
        fp_t2 = i;
    elseif flag(i) == 3
        fp_t3 = i;
    end
end


time_1 = indvar(1 : fp_t1);
time_2 = indvar(fp_t1+1 : fp_t2);
time_3 = indvar(fp_t2+1 : fp_t3);




% Define your parameters below.
k = params(1);
m = params(2);

Amp1 = params(3);
Amp2 = params(4);
Amp3 = params(5);
 
L1 = 60;
L2 = 80;
L3 = 100;
 
 
%% Section 3: Simulate Time Courses*** 
    % Here you will define your fitting equation that will be used to
    % simulate data sets for fitting. For each data set you are fitting you
    % will need a corresponding simulated data set (sim_data). This is done
    % by creating a column in sim_data for each simulated data set. Global parameters
    %are the same in all models while locals are designated with an _#(ie. A
    %for model 1 is A_1).
    
sim_data_1 = Amp1.*talbot_inversion(@(s) ((k^(L1/m))/(s*(k+s)^(L1/m))), time_1);
sim_data_2 = Amp2.*talbot_inversion(@(s) ((k^(L2/m))/(s*(k+s)^(L2/m))), time_2);
sim_data_3 = Amp3.*talbot_inversion(@(s) ((k^(L3/m))/(s*(k+s)^(L3/m))), time_3);
 
 
sim_data_compiled = [sim_data_1; sim_data_2 ; sim_data_3]; 


