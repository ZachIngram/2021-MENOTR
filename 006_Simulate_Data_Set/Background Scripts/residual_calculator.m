function [res]=residual_calculator(paropt,params,data,indvar,weights)
    %Utilizes iteration control vector from parlim to determine if each
    %parameter is being floated
for n=1:length(params(:,1))%if not, then the value doesn’t change 
    if params(n,1)==0
        params(n,3)=params(n,3);
    else %if so, then the value is allowed to change in accordance with lsqnonlin
        params(n,3)=paropt(n);
    end
end
    %Calculates residuals by subtracting the simulated data from the
    %experimental data and dividing that by the weights element wise.
res=(data-model_simulator(params(:,3),indvar))./weights;

res(~isfinite(res)) = 0; %if any element in the residual is NaN or Inf replace with 0. 