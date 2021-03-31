function [sim_data]=model_simulator(params,indvar)
%This will have to be completely edited if model is changed. The default
%layout here is an example of fitting 3 linear data sets to the equation
%y=m*x+b where each data set has the same global slope but local intercept
%values. ALL SECTIONS THAT NEED TO BE MODIFIED ARE MARKED WITH ***
 
%% Sections that require modification are:
    % Section 2 
    % Section 3
 
%% Section 1: Check Orientation 
%Corrects orientation of params matrix if needed
[r,c]=size(params);
if r>c
    params=params';
end
 
%% Section 2: Define Parameters***
% This section defines the independent variable and parameters that will be
% used in your model. Each parameter will be set equal to an element in a
% matrix called params (Example: A = params(#); where A is a fitting parameter
% and # is the position it occupies in the params matrix.)
% List all of your parameters following the example given below. These will
% be used in the model below. 
 
% If you have local parameters in your fit you will need a unique parameter
% for each data set for that parameter. See A shown below. The m parameter
% is a global parameter, but the three b values are local. 
 
    % Example:
    % m = params(1);
    % b_1 = param(2);
    % b_2 = param(3);
    % b_3 = params(4); 
 
 
% Define your independent variable below.
x = indvar;
 
% Define your parameters below.
m = params(1);          
b = params(2);        
       
 
 
 
%% Section 3: Fitting Equations***
% Here you will define your fitting equation that will be used to
% simulate data sets for fitting. For each data set you are fitting you
% will need a corresponding simulated data set (sim_data). This is done
% by creating a column in sim_data for each simulated data set. Global 
% parameters are the same in all models while locals are designated with 
% an _#(ie. A for model 1 is A_1).
    
sim_data(:,1)  = m*x+b; 
 
 
 
 
 
 
 


