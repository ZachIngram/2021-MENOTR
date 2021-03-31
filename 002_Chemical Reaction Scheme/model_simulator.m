function [Sim_Data_final]=model_simulator(params,time)
%This will have to be completely edited if model is changed. 
% ALL SECTIONS THAT NEED TO BE MODIFIED ARE MARKED WITH ***
 
% Sections that require modification are:
    % Section 2
    % Section 4 
 
 
 
 
%% Section 1: Check Inputs  
% In this section of the code, the input parameters "params" and "time" 
% are checked to make sure they're in the right orientation before moving 
% forward to the next section. 
 
[r,c] = size(params);
if r > c
    params = params';
end
[r,c] = size(time);
if r > c
    time = time';
end
 
 
tsim = [0, time]; % Creates an initial zero point
OPT = odeset('abstol',1e-12,'jacobian','off'); % ODE23tb settings
 
 
 
%% Section 2: Define Parameters, Input Initial Conditions, and Define Normalization***
% In this section of the code the user is expected to define a few inputs
% that are customized to their model. The first expectation is to define
% the parameters used in their scheme. The second expectation is to define
% the initial conditions for the model and the value for normalization.  
 
% Local Fit Setup
    k1 = params(1);                     
    k2 = params(2);   
    
    Amp = params(3);
    
    model_param = [k1 k2];
    int_conc_species = [100e-6 0];
    normalization = 100e-6;
   
    
    
% % Global Fit Setup
%     
%     k1 = params(1);                     
%     k2 = params(2); 
%     k3 = params(3); 
%  
%     L = [10e-6 100e-6 1000e-6];
%     
%     %Data Set 1    
%         param(1,:) = [k1 k2 k3];
%         int_conc(1,:) = [100e-6 L(1) 0 0];
%     %Data Set 2
%         param(2,:) = [k1 k2 k3];
%         int_conc(2,:) = [100e-6 L(2) 0 0];
%     %Data Set 3
%         param(3,:) = [k1 k2 k3];
%         int_conc(3,:) = [100e-6 L(3) 0 0];
  
 
 
 
 
%% Section 3: Simulate Time Courses 
% The simulated time courses are generated in this section. odetb is one 
% of a few built in numerical ODE solving algorithms. In this case, the 
% scheme is passed to the ODE solver along with the time points, initial 
% concentration for each species in the scheme and the model parameters. 
% The resultant output is an output matrix called sim_species_time_courses 
% where each column contains a time course for a species within the
% reaction scheme
    for i = 1:size(model_param,1)       
        
        [~,sim_species_time_courses] = ode23tb('scheme',tsim,int_conc_species(i,:),OPT,model_param(i,:));
 
%% Section 4: Define Observed Product from Reaction Scheme*** 
% In this section the product time course is identified from the matrix of
% all simulated time courses. The default in the code corresponds to the
% final species in the reaction scheme (writtedn as (:,end)). In the event
% that the product time course is composed of multiple species, the time 
% courses should be summed in this section of the code.

% For example, consider the case of A->B->C and assume the second chemical
%  species "B" in the reaction scheme corresponds to the experimentally 
% monitored signal. In this case the simu_species_time_courses is changed
% to '(:,2)' (as shown below) to correspond to the second species "B". The 
% "A" species would be '(:,1)'  and the "C" species would be '(:,3)'. 

% Ex: Product_sim_time_course(:,1) = sim_species_time_courses(:,2);


        Product_sim_time_course(:,1) = sim_species_time_courses(:,end);
 
       
%% Section 5: Normalize the Product Time Course(s) 
 
%calculate the normalized data 
Prod_num = size(Product_sim_time_course,2); % Defines the number of columns in the product time courses matrix 
        
    for j = 1:Prod_num            
        Product_sim_time_course(:,j) = Amp*Product_sim_time_course(:,j)./normalization; 
    end 
 
 
        
    
    
Sim_Data_final(1:length(time),Prod_num*(i-1)+1:Prod_num*i) = Product_sim_time_course(2:length(tsim),1:Prod_num); 
    end

