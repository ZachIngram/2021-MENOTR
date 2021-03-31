MC_sim_num = 500; 


addpath('Background Scripts')%Load Background Scripts for MENOTR
load('CompleteFit_Workspace')
[m,n] = size(data_input);
MC_sims = zeros(m,(MC_sim_num*n));


% STANDARD DEVIATION CALCULATION
% Creates a simulated data set and calculates the standard deviation of the
% simulation compared to the data
sim_data = model_simulator(param_vals_fit(1,1:end),indvar);
st_dev = sqrt(sum(sum((data_input-model_simulator(param_vals_fit,indvar)).^2))/numel(sim_data)); 


% MONTE CARLO DATA SET
    % Creates multiple sets of synthetic data points of the same size as the
    % origional data set. Each synthetic data point is picked at random from a normal
    % distribution centered about the original point. The width of the normal
    % distribution from which the new point is chosen is scaled by the standard
    % deviation calculated above. The number of sets is defined by the parameter
    % "MC_sim_num"


    
for k = 0:MC_sim_num-1
    for i = 1:m 
        for j = 1:n
            MC_sims(i,j+k*n) = st_dev.*randn(1) + sim_data(i,j);
        end 
    end
end

% WEIGHTS
    % Generates a set of weights, all of which are 1.
weights = ones(size(MC_sims)); 

param_vals = param_vals_fit;

CyMat = CyMat(end,:);
CyMatH = 1;


run_fig = 'off';
 
%calculate how much space is avaialbe in RAM
[~,sys] = memory;
Max_mem = sys.PhysicalMemory.Total;
tic
for MC_loop_num = 1:MC_sim_num
   
    %These calculations appear to require alot of RAM and reseting the
    %parallel workers has been observed to free up more memory. Here we
    %have a quick test to determine if the parallel workers should be
    %reset. If the available memory is less than 75% of the max available memory,
    %then the parallel workers will reset. 
    [mem_info] = memory;
    if mem_info.MaxPossibleArrayBytes < Max_mem *0.75 
       poolobj = gcp('nocreate');
       delete(poolobj);
   else
   end
    
    
%     clearvars -except MC_sims MC_params MC_sim_num MC_loop_num  
   
    
    MC_loop_num % displays in command window which MC simulation is being
                % fit. 
    

    data_input  = MC_sims(:,MC_loop_num);
    
    run RunAlgorithm_MC
    
    MC_params(:,MC_loop_num) = Fit_ParameterValues';
    
end

disp('Fit Complete')
disp(['Elapsed Time: ',elptim_y,' ',tunit])

%closes figures 
close all

%Calculate mean and standard deviation of the optimized parameters
%describing the MC simulated data sets
MC_mean = mean(MC_params,2);
MC_std = std(MC_params,0,2);

%Organize mean and standard deviation into output for user
MC_Mean_StDev = [MC_mean MC_std];

%Save Complete Workspace
save('Long_MC_Complete_Workspace')

clearvars -except MC_Mean_StDev

%Save parameter uncertainty values for user
save('Long_MC_Workspace')

toc