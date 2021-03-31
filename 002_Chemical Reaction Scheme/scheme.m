function ODEs = scheme(t,s,flag,k)
% ALL SECTIONS THAT NEED TO BE MODIFIED ARE MARKED WITH ***

% Sections that require modification are:
    % Section 1
    % Section 2 




%% Section 1: Define Species***
% Define the species present in your model. It is often more advantageous
% to use single letters to define species for coding purposes. 
% Example:
% A= s(1);
% B = s(2);


A = s(1);
B = s(2);


%% Section 2: Input ODE for Each Chemical Species***

% In this section write out the system of ODEs describing the scheme. 
% For the reaction A = B with forward rate constant k1 and reverse rate constant k2 the differential equations are 
% d[A]/dt = -k1*[A] + k2*[B] and
% d[B]/dt = k1*[A] - k2*[B]
% for this example the ‘ODEs' matrix is written as: 
% Example:
% ODEs = [-A*k(1) + B*k(2) 
%     A*k(1) - B*k(2)];
 
switch flag
    case ''
    
ODEs = [-A*k(1) + B*k(2)
    A*k(1) - B*k(2)];




%% Section 3 - Jacobian 

% In an effort to decrease the run time it is sometimes advantageous to
% calculate the jacobina matrix for the system of ODEs. The jacobian for
% the above reaction is shown as an example. A user supplied jacobian is 
% not necessary to execute the MENOTR fit. 
%Example: 
%  out = [-k(1) k(2)
%      k(1) -k(2)]; 


    case 'jacobian'        
 out = []; 
end
end