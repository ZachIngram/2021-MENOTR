%% GRID SEARCH ANALYSIS
    % Grid search analysis is a method of probing the constraint of a fit
    % parameter, evaluating its confidence intervals and generating a plot
    % of its error space. This is done by generating a range of values the
    % parameter will survey and floating all other parameters at each point
    % along that range to achieve a best fit(LSQNL) and corresponding chisq val.
    % From this, a plot of Chisq and parameter value can be generated, this
    % is the parameters error space. Using F statistics a F critical cut off
    % can be drawn to see how much of the error space lies below giving
    % information on the constraint of the parameter based on how narrow or
    % wide this area is. Where the F critical cut off intersects the
    % parameter error space curve gives the confidence interval for the
    % parameter.
    % In this script the user can enter which parameters are to evaluated,
    % the number of points to survey in the search and the bounds of the
    % survey range
%% Section 1: Initial Setup*** 
    clear
addpath('Background Scripts')%Load Background Scripts for MENOTR
load('CompleteFit_Workspace')
    % Initailizes which parameters the grid search will be run on. Place the
    % paramter number within the [ ].
fxparams = [2];
    % Number of points that will be surveyed above and below the best
    % fit param val.(ie 10 yeilds 21 grid search points, including the best
    % fit val)
grdpnt = 5; 
    % Scales and bounds the range that the grid search will cover. This is
    % done by scaling the best fit parameter value being evaluated by a
    % user designated percentage. As such: paramval+/-(divval*paramval) (Ex.
    % if 30% is chosen, the grid search will survey 70%-130 of the parameter
    % value) Enter the percentage as a number between 0-100
divval =50; 
 


%% Optional Parameters 
 
 
algsw =0; %Defines what NLLS solver to use. 0 is trust-region-reflective 
 % and 1 is Levenberg-Marquardt

% optional parameter value bounds that can constrain parameter values
% during parameter uncertainty calculations. 

lowerbound = 1e-5; 
upperbound = 1E10;
 
%% Section 2: Executing Grid Search 
% PARPOOL SETUP OPTIONAL
        % indexes number of available cores and sets up in parallel, this
        % is only needed when PC has more cores than matlab will use by
        % default
% delete(gcp('nocreate'))
% cluster = parcluster('local');
% numwork = cluster.NumWorkers;
% parpool(numwork)
tic %start run timer
 


 

% PREALLOCATION
ldivpnt = zeros(length(fxparams),grdpnt);
rdivpnt = ldivpnt;
ladjpnt = ldivpnt;
radjpnt = ldivpnt;
grd_fixed_param_vals = zeros(length(fxparams),(2*grdpnt+1));
grd_chisq_norm = grd_fixed_param_vals';
CyMatGrd = zeros(length(fxparams),(length(CyMat)-1));
grd_chisq = zeros((2*grdpnt+1),1);
grd_params = zeros((2*grdpnt+1),(length(CyMat)-1));
 

% Initialization of Grid Search Range and Points
div = ((divval/10)/grdpnt)/11;
for j = 1:length(fxparams)
% Generates a matrix with appropriate scaling factors to be applied to the
% best fit value for each grid point
    % Left Hand Side
    for index_fix_params = 1:grdpnt
        ldivpnt(j,index_fix_params) = (1-(div*index_fix_params));
    end
    % Right Hand Side
    for index_fix_params = 1:grdpnt
        rdivpnt(j,index_fix_params) = ((1+div*index_fix_params));
    end   
% Application of scaling factor to best fix values
    ladjpnt(j,:) = ldivpnt(j,:).*param_vals_fit(fxparams(j)); %Left hand side adjusted parameter value
    for index_fix_params = 1:grdpnt   
        if   ladjpnt(j,index_fix_params) < lowerbound 
              ladjpnt(j,index_fix_params) = lowerbound; 
       end
    end
    radjpnt(j,:) = rdivpnt(j,:).*param_vals_fit(fxparams(j)); %Right hand side adjusted parameter value
    for index_fix_params = 1:grdpnt
        if radjpnt(j,index_fix_params) > upperbound
            radjpnt(j,index_fix_params) = upperbound; 
        end
    end
    
 % Compiles all scaled parameters together
    grd_fixed_param_vals(j,:) = [ladjpnt(j,:) param_vals_fit(fxparams(j)) radjpnt(j,:)];    
 % Generation of iteration control vector for each param being
 % surveyed(current param under investigation is fixed)
    CyMatGrd(j,:) = CyMat(end,1:end-1);
    CyMatGrd(j,fxparams(j)) = 0; 
end
grd_fixed_param_vals=grd_fixed_param_vals'; % Transpose of grid search values to be used later
 
 
% Grid Search Fitting
    % Minimization and Chisq calculation of all grid search params
for index_fix_params = 1:length(fxparams) 
    disp('-------------------------------------------------------------------')
    for j = 1:2 %Division of Calculations to the right and left hand sides of the best fit value
        if j == 1
            cparam = ladjpnt(index_fix_params,:);
        else
            cparam = radjpnt(index_fix_params,:);
        end
        for kk = 1:length(cparam)
            ii = num2str(index_fix_params);
            if j == 1
            ij = 'Left';
            else
            ij = 'Right';    
            end
            ik = num2str(kk);
            disp(['Parameter:',ii,'   ','Side:',ij,'   ','Grid Search Point:',ik])
                % Defines the current parameter values for the given grid
                % seaerch point
            if kk == 1    
                param_vals_grd = param_vals_fit;
                param_vals_grd(fxparams(index_fix_params)) = cparam(kk);
            else
                param_vals_grd = grd_params(kk-1,:);
                param_vals_grd(fxparams(index_fix_params)) = cparam(kk);   
            end
            
            
            %Unique imputs for Step3 in the case of grid search analysis 
            CyMat = [CyMatGrd(index_fix_params,:) 1];
            CyMatH = 1;
            param_vals = param_vals_grd;
            g = 3;
            pop_size = 1000;
            
            run RunAlgorithm_Grid
%             
%                 % Application of iteration control vector that controls which parameters
%                 % are alound to float within the fit
%             parlim = [CyMatGrd(i,:)' params_lb' param_vals_grd' params_ub'];
%                 % LSQNL Fit
%              for r = 1:LSQNL_num    
%       
%                  [~,~,~,grd_params((((i-1)*length(cparam)*2)+((j-1)*length(cparam))+k),:),grd_chisq((((i-1)*length(cparam)*2)+((j-1)*length(cparam))+k),1),~]=lsqnl_fitting_routine(indvar,data_input,parlim,weights,algsw); 
%             parlim(:,3) = grd_params((((i-1)*length(cparam)*2)+((j-1)*length(cparam))+k),:)';
             
            grd_params((((index_fix_params-1)*length(cparam)*2)+((j-1)*length(cparam))+kk),:) = Fit_ParameterValues; 
             
           grd_chisq((((index_fix_params-1)*length(cparam)*2)+((j-1)*length(cparam))+kk),1) = Fit_ChiSq;  
            
            
             
             
             
             
             
%              end
        end  
    end    
end

%Clean up workspace before moving forward
clearvars -except grd_params grd_chisq fxparams grdpnt grd_fixed_param_vals
load('CompleteFit_Workspace')



% Chisq Organization
for index_fix_params=1:length(fxparams)
    %Rearrangment and organization of Chisq values for ease of use
grd_chisq_norm(:,index_fix_params) = [grd_chisq(((index_fix_params-1)*(2*grdpnt)+1):((2*(index_fix_params-1)+1)*grdpnt),1);best_chisq;grd_chisq((1+((index_fix_params-1)*grdpnt)+grdpnt):(((index_fix_params-1)*grdpnt)+2*grdpnt),1)];
grd_chisq_norm(:,index_fix_params) = grd_chisq_norm(:,index_fix_params)./best_chisq;
end
%F Stats Initialization
DoF = numel(data_input)-(length(CyMat)-1); % Degrees of freedom calculator
Fcrit = finv(0.68,DoF,DoF); % F critical calculator, 68% confidence(can change if needed use p)
 
 
% Confidence Interval & Plot
xql = zeros(length(fxparams),10000)';
pf = xql;
f1 = figure;
ConfIntv = zeros(2,length(fxparams));
for index_fix_params = 1:length(fxparams)
    % Interpolation of Fcrit Points and Confidence Interval
        % Interpolates a large number of data points between those that were
        % surveyed. The intersection of F crit and the error space curve
        % can be estimated from this.
xl = grd_fixed_param_vals(:,index_fix_params);
yl = grd_chisq_norm(:,index_fix_params);
xql(:,index_fix_params) = linspace(grd_fixed_param_vals(grdpnt,index_fix_params),grd_fixed_param_vals(end,index_fix_params),10000)'; % generates 10,000 x coord for interpolation points
 
%adds small scaling factor to deal with repeated values that occur as a
%consequence of bounds. 
xle = cumsum(ones(size(xl))).*xl*eps;                       % Scaled Offset For Non-Zero Elements
xle = xle + cumsum(ones(size(xl))).*(xl==0)*eps;             % Add Scaled Offset For Zero Elements
xli = xl + xle;  
 
pf(:,index_fix_params) = pchip(xli,yl,xql(:,index_fix_params)'); % interpolates the matching 10,000 y coordinate 
fc_findl = find(pf(:,index_fix_params)<Fcrit,1,'first'); % finds the left hand index for the Fcrit point(y coord)
fc_findr = find(pf(:,index_fix_params)<Fcrit,1,'last'); % finds the right hand index for the Fcrit point(y coord)
ConfIntv(:,index_fix_params) = [xql(fc_findl,index_fix_params) xql(fc_findr,index_fix_params)]'; % finds corresponding x coord for Fcrit and indexes them in a matrix
    % Plot of surveyed points, interpolation, fcrit cutoff and confidence
    % interval points
subplot(floor(2*sqrt(length(fxparams)/2)),ceil(length(fxparams)/floor(2*sqrt(length(fxparams)/2))),index_fix_params) % determines how many subplots are needed
plot(grd_fixed_param_vals(:,index_fix_params),grd_chisq_norm(:,index_fix_params),'.r','MarkerSize',20) % plots grid search points
hold on
plot(xql(:,index_fix_params),pf(:,index_fix_params),'r'); %  plots interpolation of grid search points
plot(ConfIntv(1,index_fix_params),Fcrit,'.k',ConfIntv(2,index_fix_params),Fcrit,'.k','MarkerSize',20) % plots confidence interval points
yline(Fcrit,'--k'); % plots the Fcrit cutoff
ind = num2str(fxparams(index_fix_params)); % subplot title indexes
xlabel(['Param ',ind,' Value'])
ylabel('F Calculated')
title(['Parameter: ',ind])
ylim([(Fcrit-2*(Fcrit-1)),(Fcrit+2*(Fcrit-1))]) % scales axes based on Fcrit 
end
 
% Save
    % save the elapsed time, plots and workspace
savefig(f1,'GrdSrchPlots')
close(f1)
save('CompleteGrdSrch_Workspace')
 
 
% OUTPUTS
GrdSrch_elapsedtime = toc;
if GrdSrch_elapsedtime < 60
    tunit = 's';
elseif GrdSrch_elapsedtime > 60 && GrdSrch_elapsedtime < 3600
    GrdSrch_elapsedtime = GrdSrch_elapsedtime/60;
    tunit = 'min';
else
    GrdSrch_elapsedtime = GrdSrch_elapsedtime/3600;
    tunit = 'hr';
end
elptim_y = num2str(GrdSrch_elapsedtime);
disp('Grid Search Complete')
disp(['Elapsed Time: ',elptim_y,' ',tunit])
 
 
 
% CLEAN UP
GrdSrch_Parameter_Fcalculated = grd_chisq_norm;
GrdSrch_Parameter_Values = grd_fixed_param_vals;
GrdSrch_ConfidenceInterval = ConfIntv';
GrdSrch_Fit_IndVar = xql;
GrdSrch_Fit_Fcalc = pf;
Fcritical = Fcrit;
save('GrdSrch_Workspace','GrdSrch_Parameter_Fcalculated','GrdSrch_Parameter_Values','GrdSrch_ConfidenceInterval','GrdSrch_Fit_IndVar','GrdSrch_Fit_Fcalc','Fcritical')
clear
load('GrdSrch_Workspace.mat')

