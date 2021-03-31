%% MONTE CARLO ERROR ANALYSIS
    % This script will perform a monte carlo analysis on the derived
    % fit parameters resulting in confidence probability distributions for
    % each of the derived parameters. This is achieved by taking the
    % best fit parameters for the experimental data and simulating a synthetic data
    % point for every experimental data point. This is used to calculate a
    % standard deviation for the synthetic data set that reflects the
    % quality of the fit. Some large number(MC_sim_num) of identical synthetic data
    % sets are generated using the same best fit parameters. To each of these
    % synthetic data sets is added random normally distributed noise using
    % the standard deviation as a distribution width. Each synthetic data
    % set it fit and the resulting best fit parameters are tabulated. These
    % best fits are then used to generate confidence probability
    % distributions and standard deviations for each of the parameters.
%% Section 1: Initialization***
clear
load('CompleteFit_Workspace')
addpath('Background Scripts')%Load Background Scripts for MENOTR
% Number of Monte Carlo simulation to be run(recommend 1000)
    MC_sim_num = 1000;
% Parameters to be checked for correlation. Put the parameter number in
% side the brackets []
    corrparams = [1 2];
%% Section 2: PARPOOL SETUP OPTIONAL
        % indexes number of available cores and sets up the cores in parallel, this
        % is only needed when PC has more cores than matlab will use by
        % default
% delete(gcp('nocreate'))
% cluster = parcluster('local');
% numwork = cluster.NumWorkers;
% parpool(numwork)
tic
%% Section 3: Executing MC Simulations and Calculating Parameter Uncertainty  
%PREALLOCATION
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
 
% CONTROL VECTORS
    % Creates a matrix which controls NLLS fit. First column controls which of
    % the parameters are allowed to float. This is pulled from the last
    % iteration control vector in Step2(CyMat). The second column is the lower
    % bound for each parameter. The third column is the initial guess for each
    % parameter. The fourth column is the upper bound on each parameter.
parlim = [CyMat(end,1:end-1)' params_lb' param_vals_fit' params_ub'];
 
% NLLS RUN
    % Runs the lsqnonlin script on the new synthetic data with the above
    % control vectors. This is parallelized, with each lsqnonlin routine
    % running on an individual core.
parfor k = 0:MC_sim_num-1
    k+1
[~,~,~,MC_params(:,k+1),~,~] = lsqnl_fitting_routine(indvar,MC_sims(:,n*k+1:n*k+n),parlim,weights(n*k+1:n*k+n),algsw) %#ok<PFBNS>
end
 
% CALCUALTE MC STDEV
MC_stdev = (std(MC_params,0,2));
MC_mean = (mean(MC_params,2));
 
%HISTOGRAMS
   % Generates confidence probability distributions for the derived
   % parameters that were floated during the fit. It will ignore parameters that were not
   % floated in the MC fits.
% sets up some of the plot settings
MC_his_div = MC_sim_num/10+1;
flparams = find(CyMat(end,1:end-1));
Mat_Dem = size(flparams,2);
MC_his_matx = MC_params(flparams,:);
MC_param_conf_u = (param_vals_fit'+MC_stdev);
MC_param_conf_l = (MC_mean-MC_stdev);
    % Single subplot that has the histogram for all parameters
h = figure;
for i=1:Mat_Dem
    subplot(floor(2*sqrt(Mat_Dem/2)),ceil(Mat_Dem/floor(2*sqrt(Mat_Dem/2))),i) % controls subplot size and dimensions
    histogram(MC_his_matx(i,:),MC_his_div,'Normalization','probability') % plot histogram
    hold on
    xline(MC_param_conf_u(flparams(i)),'b','LineWidth',2); % generates vertical lines tha represent 1 stdev from the mean(upper)
    xline(MC_param_conf_l(flparams(i)),'b','LineWidth',2); % (lower)
    v = param_vals_fit(flparams(i)); % indexes the mean param value
    xline(v,'r','LineWidth',2); % generates vertical line at the mean value of the parameter
    ind = num2str(flparams(i)); % indexes the param being plotted
    title(['Parameter: ',ind])
    xlabel('Best Fit Parameter')
    ylabel('Relative Probability')
    text(0.05,0.9,'- Mean','Units','normalized','FontSize',14,'Color','r')
    text(0.05,0.8,'- Std. Dev.','Units','normalized','FontSize',14,'Color','b')
end
savefig(h,'Monte Carlo Histograms') 
close(h)
    % Individual histogram is created for each parameter
for i = 1:(size(flparams,2))
    vv = flparams(i);
    h = figure;
    histogram(MC_his_matx(i,:),MC_his_div,'Normalization','probability');
    hold on
    xline(MC_param_conf_u(flparams(i)),'b','LineWidth',2);
    xline(MC_param_conf_l(flparams(i)),'b','LineWidth',2);
    v = param_vals_fit(flparams(i)); 
    plot([v v], ylim,'r','LineWidth',2)
    ind = num2str(flparams(i));
    title(['Parameter: ',ind]);
    xlabel('Best Fit Parameter')
    ylabel('Relative Probability')
    text(0.05,0.9,'- Mean','Units','normalized','FontSize',14,'Color','r')
    text(0.05,0.8,'- Std. Dev.','Units','normalized','FontSize',14,'Color','b')
    hold off
    savefig(h,sprintf('Parameter_%d_Histogram.fig',vv)) ;
    close(h)
end
 
% Parameter Correlation Plots
    % This section plots each designated parameter as a function of every
    % other designated parameter. It then generates a linear least-squares
    % fit of each along with its corresponding y-int and slope. 
for i = 1:length(corrparams)
h = figure(i);
corrparams_adj = corrparams;
corrparams_adj(i) = [];
for j = 1:length(corrparams_adj)
   subplot(floor(2*sqrt(length(corrparams)/2)),ceil(length(corrparams)/floor(2*sqrt(length(corrparams)/2))),j) % subplot size and dimension
   plot(MC_params(corrparams(i),:),MC_params(corrparams_adj(j),:),'.r','MarkerSize',5) % plot of params
   hold on
   indi = num2str(corrparams(i)); % indexes current global param
   indj = num2str(corrparams_adj(j)); % index current local param
   sgtitle(['Parameter Correlation Plots of Parameter ',indi])
   title(['Parameter: ',indj]);
   xlabel(['Parameter: ',indi])
   ylabel(['Parameter: ',indj])
   l = polyfit(MC_params(corrparams(i),:),MC_params(corrparams_adj(j),:),1); % lin lst-sqr fit of plot, generates slope and int
   xlin = 0:0.01:(MC_params(corrparams(i),1)+MC_params(corrparams(i),end)); % domain to plot fit over
   ylin = polyval(l,xlin); %generate range from lin-sqr fit
   plot(xlin,ylin,'-k') % plots fit
   slp = num2str(l(1)); % indexes slope of fit
   int = num2str(l(2)); % indexes int of fit
   text(0.6,0.95,['Slope: ' slp],'Units','normalized','FontSize',14) % displays slope and int on plot
   text(0.6,0.75,['Intercept: ' int],'Units','normalized','FontSize',14)
   ylim([ min(ylin) max(ylin)]) % controls axes limits
   xlim([min(xlin) max(xlin)])
   hold off   
end
savefig(h,sprintf('ParameterCorrelation_Parameter_%d.fig',corrparams(i))) % save each parameters plot
close(h)
end
 
% Save Fit Parameters and MC
    % Saves original best fit parameters and their corresponding stdevs
    % from the MC simulations
MC_MeanAndStDev = [param_vals_fit' MC_stdev];
save('CompleteMC_Workspace')
 
 
% OUTPUTS
MC_elapsedtime = toc;
if MC_elapsedtime < 60
    tunit = 's';
elseif MC_elapsedtime > 60 & MC_elapsedtime < 3600
    MC_elapsedtime = MC_elapsedtime/60;
    tunit = 'min';
else
    MC_elapsedtime = MC_elapsedtime/3600;
    tunit = 'hr';
end
elptim_y = num2str(MC_elapsedtime);
disp('Monte Carlo Complete')
disp(['Elapsed Time: ',elptim_y,' ',tunit])
 
 
%CLEAN UP
MC_SimualtedDataSets = zeros(size(indvar,1),1+size(MC_sims,2));
MC_SimualtedDataSets(:,1) = indvar;
MC_SimualtedDataSets(:,2:end) = MC_sims;
MC_SimulatedValues = MC_params;
save('MC_Workspace','MC_SimualtedDataSets','MC_MeanAndStDev','MC_SimulatedValues')
clear
load('MC_Workspace.mat')

