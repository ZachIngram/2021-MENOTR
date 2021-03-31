function [parlim,residual,resnorm,params,chisq,RMSD]=lsqnl_fitting_routine(indvar,data,parlim,weights,algsw)
% Function that performs the NLLS fitting portion of the algorithm, the
% actual function being minimized is located in residual_calculator. This function
% feeds in parameter values into residual_calculator and evaluates resulting
% residuals.

%% SETUP OF LSQNL:
    % nonlinear least-squares solver
OPT = optimset('lsqcurvefit'); 
    % switch that chooses which algorithm to use, this switch is set by the
    % user in step2_algorithmsetup
switch algsw
    case 0
        OPT = optimset(OPT,'Algorithm','trust-region-reflective');
    case 1
        OPT = optimset(OPT,'Algorithm','levenberg-marquardt');
end
    % When possible the algorithm will approximate gradients in parallel
OPT = optimset(OPT,'UseParallel',false);
    % Termination tolerance of function
OPT = optimset(OPT, 'Tolfun', 1e-10);
    % Termination tolerance on the current point
OPT = optimset(OPT, 'TolX', 1e-10);
    % Maximum number of iterations or function evaluations, whichever comes first
% OPT = optimset(OPT,'MaxIter',50,'MaxFunEvals',100);
    % Results of output is not echoed
OPT = optimset(OPT,'display','off');

%% LSQNONLIN FUNCTION:

switch algsw
    case 0
        [x,resnorm,wres,~,~,~,~]=lsqnonlin('residual_calculator',parlim(:,3),parlim(:,2),parlim(:,4),OPT,parlim,data,indvar,weights);    
    case 1
        [x,resnorm,wres,~,~,~,~]=lsqnonlin('residual_calculator',parlim(:,3),[],[],OPT,parlim,data,indvar,weights);
end
%Redefines current parameters with output from the NLLS algorithm    
parlim(:,3)=x;
    %Calculates and defines all of the outputs for the overall lsqnl_fitting_routine
    %function
residual=wres;
RMSD=sqrt(sum(sum(residual.^2))/numel(data));
chisq=sum(sum(wres.^2));
params=parlim(:,3);
