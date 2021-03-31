%% ALGORITHM SETUP AND INITIALIZATION
% This script will be the most heavily modified section of the code. It
% initializes the fit settings and sets up the data and initial guesses. ALL
% SECTIONS THAT NEED TO BE MODIFIED ARE MARKED WITH ***
 
%% Sections that require modification are:
    % Section 1
    % Section 2 
    
 
 
 
 
%% Section 1: PARAMETER GUESSES***
    % Defines the initial guesses that the algorithm will use as a starting
    % point. The GA portion of the algorithm will help reduce how sensitive
    % the fit is to the initial guesses; however, these guesses need to be
    % reasonable(preferably within an order of magnitude). Each parameter
    % is also assigned an upper and lower bound. Note: Upper and Lower
    % bounds have to be different number
    % Param 1
 
k1 = 10;                  % Initial guess
k1_ub = 100;               % Upper bound
k1_lb = 0.001;               % Lower bound
 
k2 = 10;                 
k2_ub = 100;                
k2_lb = 0.001; 
 
Amp = 1; 
Amp_ub = 1;
Amp_lb = 0.01;
 
   
 
    % Concatenates params and their bounds into a single matrix ***
    % keep the format of [x x x x x x...];
param_vals = [k1 k2 Amp];
params_ub = [k1_ub k2_ub Amp_ub];
params_lb = [k1_lb k2_lb Amp_lb];
clear k1 k2 Amp k1_ub k2_ub Amp_ub k1_lb k2_lb Amp_lb
 
%% Section 2: FIT SETTINGS***
% Iteration Control Vectors:
    %Controls how many cycles the algorithm will perform. Mute unneeded
    %cycles and create more by copying and editing the index number to
    %match the cycle number(Cycle #** --> kt_c(**)=1;). Within each cycle
    %set the number of generations by editing the value of g and control
    %which parameters will be floated in each cycle(1=float 0=don’t float)
    %Set LSQNL# equal to 1 to turn on, set equal to 0 to turn off 
%RECOMMENDATION: 
%The default iteration control vector should consist of a
%first cycle that floats all of the fit parameters for 10 generations with
%LSQNL turned off, with a second cycle that floats all the parameters for
%10 generations with LSQNL turned on. This works for most cases.
 
 
                             % Cycle #1
 
g(1)=5;     % g(1) is the first element in the generations matrix which  
            % corresponds to the number of generations to be performed 
            % in the first cycle, 5 indicates how many generations to 
            % perform within the first cycle. 
            
k1_c(1) = 1;  % m_c(1) : m is the parameter name, _c indicates the cycles 
              % matrix, (1) indicates the first cycle, = 1 indicates this
              % is a floating parameter. A floating parameter whose value
              % is allowed to change within the model. A parameter which is
              % not allowed to float is not allowed to change and is
              % instead a fixed parameter value. 
k2_c(1) = 1;
Amp_c(1) = 1;
 
 
NLLS(1)=0; % This defines if the first cycle should perform NLLS (= 1) 
            % or not (= 0). It is generally suggested that NLLS not 
            % be performed on the first cycle. 
 
            
                    
            
                             % Cycle #2
                             
g(2)=5;     % g(2) is the second element in the generations matrix which 
            % corresponds to the number of generations to be performed in 
            % the second cycle in the second ...
            % cycle, 3 indicates how many generations to perform within
            % the second cycle.  
 
k1_c(2) = 1; 
k2_c(2) = 1;
Amp_c(2) = 1;
 
 
NLLS(2)=1;
 
    % Concatenates iteration control vector***
    % format as such [x;x;x;x;x]; (need ; between each term)
    
CyMat=[k1_c; k2_c; Amp_c; NLLS]';
CyMatH=size(CyMat,1);
clear k1_c k2_c Amp_c NLLS
 
%% Section 3: RUN FIGURE SETTINGS
    % Run_fig
        % This controls whether a run figure is generated throughout the
        % fitting process. This figure shows general fit information, the
        % current fit of the data, and the change in chisq. Set to 'on' or 'off'
    run_fig = 'on';
    % Axspace
        % Choose whether the data displayed in the run figure is on a
        % linear or log scale. Linear 'lin' or Log 'log'
    axspace = 'log';

 
%% Section 4: DATA INDEXING
 
    %Partitions data into the independent variable, dependent
    %variable and the StDev on the dependent variable.
indvar = InsertData(:,1); % pulls out indvar
data_input  = InsertData(:,2:2:end); %pulls out data
weights = InsertData(:,3:2:end); %pulls out weights
 
for i=1:size(weights,1) %weights should not be equal to zero. This checks weights matrix. 
    for j=1:size(weights,2)
        if weights(i,j)<0.00001
            weights(i,j)=0.00001;
        end
    end 
end
 
 
 
%% Section 5: OVERALL SETTINGS
    % Stopcrit
        % Defines the stopping criteria. The algorithm will conclude it has
        % converged and stop when (1) the mean, absolute coefficient of variance(cv) of the
        % parameters is less than the stopcrit_params (2) the percent
        % difference between the chisq at the beginning of period and the
        % best chisq at the end of the period is less than the
        % stopcrit_chisq (3) the fit has run stopcrit_period periods(this
        % last stopping criteria is to ensure the fit doesn’t get stuck and
        % go on indefinitely) 
    stopcrit_param = 0.01; % cutoff for parameter cv (0.01 is recommended default)
    stopcrit_chisq = 1; % cutoff for percent difference of chisq(1% is recommended default)
    stopcrit_period = 10; % cutoff for maximum number of periods (10 is recommend default)
    % ORnum
        % Defines the number of optimization routines that each period will run through(3 is the recommended default).
     ORnum = 3;
    % Pop_size
        % Indicates number of individuals within each generation(for 1-5
        % parameters 20 is the recommended default, 6-10=200, 11-15=1000,
        % >15=2000)
    pop_size = 200; 

    
    
    % GA SETTINGS:
    % Crossover
        % Fraction of individuals that are generated from crossover in each
        % generation(rest of the fraction is generated from mutation)(0.85 is the recommended default)
    crossover = 0.85;
    % Stdevs
        % Defines the normal distribution by which random mutations are
        % generated(the is fractional, 1=100% of value of the parameter)(1 is the recommended default)
    stdevs = ones(size(param_vals));
    % Best
        % Number of individuals used to make the best population.(must be
        % greater than nlls_size) (5 is the recommended default)
    best = 5;
    % Mut_probability
        % This number defines the probability that a parameter value within
        % an individual will be mutated.  When individuals are being
        % produced by mutation, the probability that the value of parameter
        % (n) of a specific individual will be produced by mutation is
        % mut_probability. The corollary is that the probability of the
        % value of parameter (n) of a specific individual running through
        % the mutation process will reach the next generation un-altered is
        % 1-mut_probability. (0.5 is the recommended default)
    mut_probability = 0.5;
    % Elite
        % This is the number of individuals that advance to the next
        % generation unaltered.  The individuals are chosen by ranking the
        % Chisq values.  The individuals that advance are the top (n) where
        % the value of "n" is defined by the elite field.  The way the
        % algorithm is written, elite must be less than or equal to best. (3 is the recommended default)
    elite = 3;
    % Decay
        % This controls the rate of change of the standard deviations used
        % during mutation.  The standard deviation changes according to the
        % reciprocal of the generation number according to
        % stdev(generationN)=stdev(generation1)*abs((1-decay*(generationN/generationT)))
        % where generation "N" refers to the generation number and
        % generation "T" refers to the total number of generations.
        % Depending on the value of decay the standard deviation changes in
        % different ways.  For decay of zero, standard deviation is
        % constant.  For decay of one, standard deviation shrinks from the
        % initial standard deviation to zero at a rate of 1/generation. For
        % a decay between 1 and generations, the standard deviation shrinks
        % at a rate dictated by decay until decay equals the generation
        % number then grows.  For decay values larger than generations the
        % standard deviation grows throughout the run. (1 is the recommended default)
    decay =1;

    
    % NLLS SETTINGS:
    % AlgSw
        % Sets which algorithm the lsqnl will utilize. Default setting is
        % 0 which uses 'trust-region-reflective'. The other option is 1 for
        % 'levenberg-marquardt'. For more info type "doc Choosing the
        % Algorithm" into the command window and look under lsqnonlin
    algsw = 0;
    % Nlls_size
        % Defines the number of individuals that will be refined via NLLS analysis
        % in each generation where NLLS is being performed.(must be lower
        % than best) (3 is the recommended default)
    nlls_size = 3;

    %% SAVE SETTINGS
    
save('Start_Workspace.mat')

