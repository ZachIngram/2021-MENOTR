%% RUN ALGORITHM
% This script should simply be ran after running Step2_AlgorithmSetup. 
 
 
%% Sections that require modification are:
    % None. A user should not edit this script. 
    
 
 
%% PARPOOL SETUP
        % indexes number of available cores and sets up in parallel, this
        % is only needed when PC has more cores than matlab will use by
        % default
%     delete(gcp('nocreate'))
%     cluster = parcluster('local');
%     numwork = cluster.NumWorkers;
%     parpool(numwork);
% addpath('Background Scripts')%Load Background Scripts for MENOTR
tic % start timer
%% LOAD ALGORITHM SETUP
% diary MENOTR_Fit_Command_Window
% load('Start_Workspace')
%% PREALLOCATE AND INITIALIZATION
    %Preallocate
param_vals_or = ones(ORnum,length(param_vals));
chisq_or = ones(ORnum,1);
best_chisq_or = zeros(stopcrit_period,1);
delta_chisq = zeros(stopcrit_period,1);
param_cv = zeros(stopcrit_period,1);
saved_period_params = zeros(stopcrit_period,size(param_vals,2));
    %Initialize
init_or_params = param_vals; %sets initial optimization routine parameters as the intial guess
stopflag = 1; % Initialize stop criteria
Period = 1;   % Period counting index
%% PERIOD AND CYCLE CONTROL        
    % Period loop
        % This first loop controls the highest level of iteration, the
        % period. Within each period, the routine will successively run
        % through some number of optimization routines(OR) defined by
        % ORnum. The first OR will start with the initial parameter guesses
        % and work through the designated number of cycles(lower level of
        % iteration). The routine will save the results of the first OR and
        % initialize the second OR using the same initial guesses as the
        % first. This continues until the first period completes when all
        % ORs are run. At the end of each period the algorithm will
        % calculate the mean coefficient of variance of the parameters
        % across all ORs and the percent difference between the starting
        % chisq of the period and the best final chisq of the period. If
        % these are lower than the stop criteria or the maximum number of
        % periods is reached then the algorithm will stop, if not the
        % algorithm will continue to the next period. The fit parameters
        % from each OR are ranked and the best is chosen to be passed to
        % the next period or saved as an output if the algorithm stops. The
        % next period will use those best fit parameters as the initial
        % guesses for all of its ORs. This repeats until these criteria are
        % met and it saves the best fit parameters from the last period as
        % the final fit parameters.
while stopflag == 1
    param_vals = init_or_params;
        % OR loop
            % This second loop controls the next level of iteration, the
            % optimization routine(OR). Within a period some number of ORs(
            % this is controlled with a user variable OR_num) are run in
            % successive order independent of one another. Each OR will
            % proceed through the two lowest tier fitting iterations,
            % cycles and generations producing an estimate of the fit
            % parameters. These ORs are run in this manner as a way to
            % access the degree of convergence the algorithm has reached
            % and increase the routines search space. Periods that contain
            % ORs with different estimates of parameters are not fully
            % converged while a period that results in all ORs producing
            % the same parameter outcomes has converged.
    for i = 1:ORnum
            % Cycle loop.
                % This third loop controls the second to last level of
                % iteration, the cycle. Within each cycle a selected group
                % of parameters are allowed to float within the fit while
                % the rest are held constant, this is controlled by the
                % user. This allows more constrained and less correlated
                % parameters to float to decent values without the less
                % constrained and more correlated parameters diverging
                % the fit. The user also has the ability to tun the nlls
                % portion of the code on and off for each cycle. These two
                % abilities allows for the user to set up a scenario where
                % a few initiation cycles can run using only the genetic
                % algorithm and generate at set of good guesses to be used
                % in the last few cycles by the nlls solver to reach an
                % optimum value. This loop has two cases, (1) CyMatH>1, meaning
                % there is more than one cycle. Here it first runs through
                % all but the last cycle not saving the running parameters,
                % then runs the last cycle and saves the needed parameters.
                % (2) CyMatH=1, there is only one cycle to run. It will
                % run and save all the needed parameters.
        if CyMatH > 1 
            for k = 1:CyMatH-1
            zy = 0;% Counting index for user display 
                [param_vals(k+1,:)] ...
                    = ga_fitting_routine_Long(param_vals(k,:),g(k),pop_size,data_input,weights,...
                crossover,stdevs,best,mut_probability,elite,decay,nlls_size,indvar,CyMat(k,:),...
            zy,g(k),params_ub,params_lb,k,CyMatH,Period,i,algsw,axspace,run_fig);
            end
            zy = 1; k = 1;% Counting index for user display                
                [param_vals(CyMatH+1,:),best_chisq,chi_values,RMSD,indvar1_sim,simvals]...
                     = ga_fitting_routine_Long(param_vals(CyMatH,:),g(CyMatH),pop_size,data_input,weights,...
                crossover,stdevs,best,mut_probability,elite,decay,nlls_size,indvar,CyMat(CyMatH,:),...
                zy,g(CyMatH),params_ub,params_lb,k,CyMatH,Period,i,algsw,axspace,run_fig);
        else
            zy = 1; k = 1;% Counting index for user display              
                [param_vals(CyMatH+1,:),best_chisq,chi_values,RMSD,indvar1_sim,simvals]...
                    = ga_fitting_routine_Long(param_vals(CyMatH,:),g(CyMatH),pop_size,data_input,weights,...
                crossover,stdevs,best,mut_probability,elite,decay,nlls_size,indvar,CyMat(CyMatH,:),...
                zy,g(CyMatH),params_ub,params_lb,k,CyMatH,Period,i,algsw,axspace,run_fig);
            
            % End of Cycle loop
        end 
        % Creates matrix of best fit for each OR in a period
    param_vals_or(i,:) = param_vals(end,:);
        % Resets param_vals back to the initial guess
    param_vals = init_or_params;
% End of OR loop
    end
        % Determines which of the ORs generated the best fit parameters and
        % indexes them
    for n = 1:ORnum
        chisq_or(n,1) = sum(sum(((data_input-model_simulator(param_vals_or(n,:),indvar))./weights).^2));
    end
    [~,J] = min(chisq_or); % sorts the chisq and indexes there position realative to their corresponding parameters
    best_chisq_or(Period,1) = chisq_or(J,1); % sets the best chisq equal to the 1st in the sorted chisq(lowest value) 
    % command window output at end of period
    bestchisq_y = num2str(chisq_or(J,1));
    per_y = num2str(Period);
    disp('-------------------------------------------------------------------')
    disp(['Period ',per_y ,' Chisq Results:'])
    disp(chisq_or)
    disp(['Best Chisq:  ',bestchisq_y])
    disp('-------------------------------------------------------------------')
    %% Stopping criteria:
        % There are two criteria that have to be meet before the algorithm
        % will stop. The first is based on the idea that each optimization 
        % routine should result in the same chi-squared. The coefficient of 
        % variation is calculated for the pooled chi-squared values and then 
        % compared to a user supplied tolerance. 
        % The second compares the mean, absolute coefficient of
        % variance of the best parameter values from each OR. When both
        % of these criteria fall below the user decided cutoff the fit will
        % stop. There is also a trigger for stopping the fit if it runs
        % over the maximum number of periods.
 
        % calculates the coefficient of variation (normalized standard
        % deviation) between the differen optimization routine chi-squared 
        % values. The coefficient of variation must be below a threshold to
        % satisfy 1 of two stopping criteria. 
        
    CV_chisq(Period,1) = std(chisq_or)/mean(chisq_or);
    
    CV_chisq(isnan(CV_chisq))=0;
    
    %calculates the mean, absolute coefficient of variation for the
        %best parameters from each OR
    param_cv_std = std(param_vals_or,0,1); % standard deviation of each parameter
    param_cv_std(~param_cv_std) = 1E-50; % replaces all zeros with 1E-30 so there is no Inf or NaN
    param_cv_mean = mean(param_vals_or,1); % mean of each parameter
    param_cv_mean(~param_cv_mean) = 1E-20; % replaces all zeros with 1E-20 so there is no Inf or NaN
    param_cv(Period,1) = mean(abs(param_cv_std./param_cv_mean)); % mean, absolute cv
        %command window output
    pcv_y = num2str(param_cv(Period));
    dcs_y = num2str(CV_chisq(Period)*100);    
    totper_y = num2str(stopcrit_period);
    disp('  ')
    disp(['Relative Standard Deviation of Period Parameters: ',pcv_y])
    disp(['Relative Standard Deviation of Chi-Squared Between ORs: ',dcs_y,'%'])
    disp(['Period ',per_y,' out of the max ',totper_y])
    disp('  ')
        % compares stopping criteria to the fit results
    if param_cv(Period) <= (stopcrit_param) & CV_chisq(Period) <= (stopcrit_chisq/100) | Period == stopcrit_period %#ok<OR2,AND2>
        stopflag = 0; %algorithm stops
    else
        stopflag = 1; %algorithm continues
    end
%     save('Period_Workspace') % save the results of the period
    if stopflag == 1
        Period = Period+1; % period counter
    end 
        % Passes the best fit params to the next period as its initial guesses
    init_or_params(1,:) = param_vals_or(J,:); 
        % Saves the best fit parameter values for each period
    saved_period_params(Period,:) = param_vals_or(J,:); 
%     save('param_vals_period','saved_period_params'); 
    % End of period loop. 
end
%% TRIMING AND CLEARING
    %Trims preallocated matrices to actual size
best_chisq_or = best_chisq_or(Period,1);
CV_chisq = CV_chisq(Period,1);
param_cv = param_cv(Period,1);
saved_period_params = saved_period_params(Period,size(param_vals,2));
%% OUTPUTS
    % elapsed time of fit
Fit_elapsedtime = toc;
if Fit_elapsedtime < 60
    tunit = 's';
elseif Fit_elapsedtime > 60 & Fit_elapsedtime < 3600
    Fit_elapsedtime = Fit_elapsedtime/60;
    tunit = 'min';
else
    Fit_elapsedtime = Fit_elapsedtime/3600;
    tunit = 'hr';
end
elptim_y = num2str(Fit_elapsedtime);
    % param vals
param_vals_fit = param_vals_or(J,:);
    %saves entire work space
% save('CompleteFit_Workspace')
disp('Fit Complete')
disp(['Elapsed Time: ',elptim_y,' ',tunit])
    % workspace clean up and plot
Sim_Fits = model_simulator(param_vals_fit,indvar);
Sim_IndVar = indvar;
Sim_Data = zeros(size(Sim_IndVar,1),(size(Sim_IndVar,2)+size(Sim_Fits,2)));
Sim_Data(:,1) = Sim_IndVar;
Sim_Data(:,2:end) = Sim_Fits;
Residual = data_input-model_simulator(param_vals_fit,indvar);
Fit_Residual = zeros(size(indvar,1),(size(indvar,2)+size(Residual,2)));
Fit_Residual(:,1) = indvar;
Fit_Residual(:,2:end) = Residual; 
Fit_Variance = sum(sum(Residual.^2))/numel(data_input);
Fit_RMSD = sqrt(Fit_Variance);
Fit_ChiSq = best_chisq;
Data = InsertData;
Fit_ParameterValues = param_vals_fit;
% FINAL PLOT
switch axspace
    case 'lin'
final_plot = figure(1);
subplot(2,1,1)
plot(indvar,data_input,'o')
hold on
plot(Sim_IndVar,Sim_Fits,'k')
title('Fits of Data')
xlabel('Independent Variable')
ylabel('Dependent Variable')
hold off
subplot(2,1,2)
plot(indvar,Residual)
title('Fit Residuals')
xlabel('Independent Variable')
ylabel('Residuals')
    case 'log'
final_plot = figure(1);
subplot(2,1,1)
semilogx(indvar,data_input,'o')
hold on
semilogx(Sim_IndVar,Sim_Fits,'k')
title('Fits of Data')
xlabel('Independent Variable')
ylabel('Depdendent Variable')
hold off
subplot(2,1,2)
semilogx(indvar,Residual)
title('Fit Residuals')
xlabel('Independent Variable')
ylabel('Residuals')
end
% SAVE
% saveas(final_plot,'Final_Plot.tif');
% save('Fit_Workspace','Data','Fit_ChiSq','Fit_ParameterValues','Fit_Residual','Fit_RMSD','Fit_Variance','Sim_Data')
% clear
% load('Fit_Workspace.mat')
% delete('Period_Workspace.mat','param_vals_period.mat','CurrentParamVals.mat')
diary off

