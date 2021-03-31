%This script is used to simulate a data set with user defined parameter
%values. 

% Sections that required editing for new simulated data set

% Section 2
% Section 3
% Section 4 

%% Section 1: Add Bakground Scripts to Working Directory 

 addpath('Background Scripts')%Load Background Scripts for MENOTR
 
%% Section 2: Define the independent variable*** 
indvar_start = 1;
indvar_end = 100;
points = 100;
 
 
%% Section 3: Identify if the independent variable should be linear or log spaced*** 
 
%%% Linearly Space Data points 
% indvar = linspace(indvar_start,indvar_end,points);
 
%%% Logarithmically Space Data points 
indvar = logspace(log10(indvar_start),log10(indvar_end),points);
 
%% Section 4: Define Parameter Values*** 
 
 
%%parameter info (this must be set up to match your model_simulator script,
%%i.e. the parameters have to be in the same order as how they are shown in
%%the model_simulator script)
 
A1 = 1;
k1 = 1;
 
 
 
 
%The params vector shown below must be edited to account for all your
%parameters, and the parameters must be in the correct order as shown in
%your model_simulator script. 
params =[A1 k1];
 
 
 
%% Section 5: Create simulated data set

%Portion of code that will create your simulation  
 
 
sim_data = model_simulator(params,indvar);
 
Sim_Data_Set = [indvar',sim_data];
 
clearvars -except Sim_Data_Set
 
%The workspace contains a matrix called "Sim_Data_Set". The first column 
%in this matrix is the independent variable while the remaining columns
%correspond to the simulated data set. 
 


