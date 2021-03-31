function [sim_data,Xf]=model_simulator(params,Xt)
%% Sections that require modification are:
    % Section 2
 
 
 
%% Section 1: Check Inputs  
% In this section of the code, the input parameter "params"  
% is checked to make sure they're in the right orientation before moving 
% forward to the next section.
%Corrects orientation of params matrix if needed
[r,c]=size(params);
if r>c
    params=params';
end
 
%% Section 2: Define Xbar, Mt values, and update code
% Defines values from params matrix into their corresponding variable 
% which will be used in the model below.
 
 
Xbar = 'Amp*(K1*Xf+2*K1*K2*Xf^2)/(1+K1*Xf+K1*K2*Xf^2)'; 
 
Mt(1) = 1E-6;

 
% Define Parmaeters
Amp = params(1);
K1 = params(2);
K2 = params(3);
 
% This code converts the parameter value from a number to a
% string. Edit the code below to match your parameters. 
Amp_s = num2str(Amp);
K1_s = num2str(K1);
K2_s = num2str(K2);
 
% This code inputs the parameter values into the Xbar equation. Edit 
% the code to match your parameters. 
Xbar = strrep(Xbar,'Amp',Amp_s);
Xbar = strrep(Xbar,'K1',K1_s);
Xbar = strrep(Xbar,'K2',K2_s);
 
 
 
 
 
 
%% %% Section 3: 
    %There is a model for each of the different data sets. Global parameters
    %are the same in all models while locals are designated with an _#(ie. A
    %for model 1 is A_1).
    
    NumMt = size(Mt,2);
    
    sim_data = ones(size(Xt,1),NumMt); %Pre-allocation of matrix  
    
    
    for k=1:NumMt
        % This section feeds Xt(1) in as the initial starting point for the
        % first run of the root finder and after that it will feed back in
        % the previous Xf as the new initial starting point. The hope is
        % this will make it more likely that the root finder is outputting
        % the wanted root.
        for i=1:size(Xt,1)
            [sim_data(i,k),Xf]=RootFinder(Xt(i),Mt(k),Xbar);   
       
        
        end
    end


