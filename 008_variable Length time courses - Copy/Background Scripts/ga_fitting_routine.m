function [param_vals,best_chisq,chi_values,RMSD,indvar1_sim,simvals,flag_sim] = ga_fitting_routine(params,generations, ...
    pop_size,data,weights,crossover,stdevs,best,mut_probability,elite,decay,nlls_size,indvar,CyMat,zy,gzy, ...
    params_ub,params_lb,kcy,CyMatH,Period,ORnum,algsw,axspace,run_fig,flag)
 
 
 
%% GENETIC ALGORITHM AND NLLS FITTING ROUTINE:
    %This script runs the lowest level of iteration, the generation. It
    %sets up the initial population at the beginning of each generation.
    %Applies the boundary conditions and iteration control vectors. If
    %called for it will initiate the nlls routine on a portion of the best
    %population. It then sets up the genetic operators and applies them to
    %the population. It also has the capability to generate a user display
    %that outputs the data with the current fit, history of the chisq, the
    %current and previous fit parameters, as well as a count of the
    %generation, cycle, and period. This display is turned off by default
    %because it isn’t general and wont work for all model. It can be edited
    %if a display is desired.    
%% INITIAL SETUP FOR A GENERATION:
    % Defines the range of the simulated data set
    
for i = 1: length(flag)
    if flag(i) == 1
        fp_t1 = i;
    elseif flag(i) == 2
        fp_t2 = i;
    elseif flag(i) == 3
        fp_t3 = i;
    
    end
      
end

time_1 = indvar(1 : fp_t1);
time_2 = indvar(fp_t1+1 : fp_t2);
time_3 = indvar(fp_t2+1 : fp_t3); 


    
    
    
    switch axspace
        case 'lin' % linear spacing
            indvar1_sim_1 = linspace(min(time_1),max(time_1),1e3)';
            indvar1_sim_2 = linspace(min(time_2),max(time_2),1e3)';
            indvar1_sim_3 = linspace(min(time_3),max(time_3),1e3)';        
        
        case 'log' % log spacing
            indvar1_sim_1 = logspace(log10(min(time_1)),log10(max(time_1)),1e3)';
            indvar1_sim_2 = logspace(log10(min(time_2)),log10(max(time_2)),1e3)';
            indvar1_sim_3 = logspace(log10(min(time_3)),log10(max(time_3)),1e3)';
    end
    
    indvar1_sim = [indvar1_sim_1; indvar1_sim_2; indvar1_sim_3];
    flag_sim = [1*ones(1000,1); 2*ones(1000,1); 3*ones(1000,1)]; 
    
    
    % Preallocated matrices for the population matrix, the "best
    % population" matrix, Chi-square matrix, saved parameters, chisq
    % values, and NLLS_pop
population_mat = zeros(pop_size,length(params));
best_pop = zeros(best,length(params));
chisq = population_mat(1:end,1);
NLLS_pop = zeros(nlls_size,length(params));
    %Prellocates and defines running variables that are used within the
    %fitting routine
running_chi = zeros(generations,1);
    % These lines generate sub-matricies for the mutation and crossover
    % population and the elite and nlls populations.  This will ensure that
    % the mutation and crossover fractions are preserved.
if nlls_size>0
    mut_cross_size = pop_size-(elite+nlls_size+1);
    elite_and_nlls_mat = zeros(elite+nlls_size+1,length(params));
else
    mut_cross_size = pop_size-elite;
    elite_and_nlls_mat = zeros(elite,length(params));    
end
mut_cross_mat = zeros(mut_cross_size,length(params));
chi_values = ones(pop_size,generations);
    % This is the initial production of a randomized population.  All
    % individuals are produced by a mutation event in which a random value
    % drawn from a gaussian distribution (mean=initial guess parameter,
    % stdev=initial guess parameter value) is added to every parameter of
    % every individual. After this saturating mutation event, a single
    % individual from the population matrix is replaced by the initial
    % guess individual.
for i = 1:length(params)
    parfor j = 1:pop_size
        params_log = log10(params)
        population_mat(j,i) = normrnd(params_log(i),1);
        population_mat(j,i) = 10.^(population_mat(j,i));
        population_mat(j,i) = real(population_mat(j,i));
    end
end
population_mat(1,1:end) = params;
%% MAIN GENERATION LOOP
    % This is the beginning of the main loop that will be executed the
    % number of times dictated by generations.  At every execution of the
    % loop the standard deviation shrinks according to
    % stdev(generationN)=stdev(generation1)*(1-decay*(generationN/generationT))
    % where generation refers to the generation number and generation
    % refers to the total number of generations.  At the last generation
    % stdev=0.
for g = 1:generations
       if g>1
       params = param_vals; 
        for i = 1:length(params)
            parfor j = 1:pop_size
                    params_log = log10(params)
                    population_mat(j,i) = normrnd(params_log(i),1);
                    population_mat(j,i) = 10.^(population_mat(j,i));
                    population_mat(j,i) = real(population_mat(j,i));
            end
        end
        population_mat(1,1:length(params)) = params; 
        end  
    
    shrink_vals = stdevs.*abs((1-(decay*(g/generations))));
        % command window display
    g_y = num2str(g);
    or_y = num2str(ORnum);
    s_y = num2str(Period);
    if zy == 0
        c_y = num2str(kcy);
        else
        c_y = num2str(CyMatH);
    end    
disp(['Period:',s_y,'   ','OptRoutine:',or_y,'   ','Cycle:',c_y,'   ','Generation:',g_y])    
    %% BOUNDS ON PARAMETERS:
        % The following loop imposes the bounds on each parameter value.
    for k = 1:pop_size
        for i = 1:(size(CyMat,2)-1)
             if population_mat(k,i) >= params_ub(i)
                 population_mat(k,i) = params_ub(i);
             elseif population_mat(k,i) <= params_lb(i)
                 population_mat(k,i) = params_lb(i);
             else
                 population_mat(k,i) = population_mat(k,i);
             end
        end
    end                     
    %% ITERATION CONTROL:
        % Utilizes iteration control vectors to set which parameters will be
        % fixed or floated in the generations for a given cycle
    for k = 1:pop_size
        for i = 1:(size(CyMat,2)-1)
             if CyMat(i) == 0
                 population_mat(k,i) = params(i);
             else
                 population_mat(k,i) = population_mat(k,i);
             end
        end
    end           
    %% CHISQUARED CALCULATION:         
        % Chisq values are calculated in parallel in the following parfor loop.
    parfor q = 1:pop_size
        chisq(q,1) = sum(sum(((data-model_simulator(population_mat(q,1:end)',indvar, flag))./weights).^2));  %#ok<PFBNS>
    end
        % pop_mat is generated for use downstream.  pop_mat is equal to the
        % population matrix from the generation that was just used to calculate
        % Chisq values.  Chisq values are sorted (ascending order) and two vectors are generated
        % based on this sorting process.  One vector (Y) is the sorted Chisq
        % values, the other vector (index) is a vector listing the positions of
        % the original Chisq vector. Index reports on the coordinates of the
        % sorted Chisq values
    pop_mat = population_mat;
    [~,index] = sort(chisq);
    %% BEST POPULATION:
        % This loop produces the "best" population.  The number of individuals in
        % the best population is defined in the input in the "elite" field.
    for m = 1:best
        best_pop(m,:)=population_mat(index(m),:);
    end
    %% NLLS Routine
    if nlls_size > 0
        % Creates a matrix which controls NLLS fit. First column controls which of
        % the parameters are allowed to float. This is pulled from the 
        % iteration control vector.  The second column is the lower
        % bound for each parameter. The third column is the initial guess for each
        % parameter. The fourth column is the upper bound on each paramter. 
        parlim = [CyMat(end,1:end-1)' params_lb' ones(length(params),1) params_ub'];
        % This next block generates the "NLLS Pop".  In the parfor loop the top
        % of the best population is further optimized by NLLS.  The number of
        % individuals to be optimized is specified by the "nlls_size" field in
        % the input. 
        if CyMat(length(CyMat)) == 1   
            parfor r = 1:nlls_size
                loop_lim = parlim;
                pars = best_pop(r,:)';            
                loop_lim(1:length(params),3) = pars(1:length(params));                        
                [fit_parlim,~,~,~] = lsqnl_fitting_routine(indvar,data,loop_lim,weights,algsw,flag);
                NLLS_pop(r,:) = fit_parlim(:,3);                
            end
        end
    end
    %% GENETIC ALGORITHM OPERATORS:%%% 
    %CROSSOVER AND MUTATION:
        % This loop generates the half of the crossover fraction that is
        % generated from the best population.
    for n = 1:round(mut_cross_size*crossover*0.5)
        for j = 1:length(params)
            mut_cross_mat(n,j) = best_pop(randi(best),j); % Uniform multi-parent approach, each gene of the new individual is take from a random individual in the best_pop
        end
    end
        % This loop generates the half of the crossover fraction that is
        % generated from the general population.
    for n = round(mut_cross_size*crossover*0.5)+1:round(mut_cross_size*crossover)
        for j = 1:length(params)
                mut_cross_mat(n,j) = pop_mat(randi(pop_size),j); % Uniform multi-parent approach
        end
    end
        % This loop generates the half of the mutation fraction that is
        % generated from the best population.
    for r = round(mut_cross_size*crossover)+1:round(length(mut_cross_mat(:,1))*0.75)
            rand_coordinate=randi(best);
        for j = 1:length(params)
            rand_num = rand();
            if rand_num <= mut_probability % selects which genes to mutate
                mut_cross_mat(r,j) = normrnd(best_pop(ran_coordinate,j),abs(best_pop(rand_coordinate,j))*shrink_vals(j)); % mutates gene based on original value and scales by shrink_vals
            else
                mut_cross_mat(r,j) = best_pop(rand_coordinate,j);
            end
        end
    end
        % This loop generates the half of the mutation fraction that is generated
        % from the general population.
    for r = round(length(mut_cross_mat(1:end,1))*0.75)+1:mut_cross_size
        rand_coordinate = randi(pop_size);
        for j = 1:length(params)
            rand_num = rand();
            if rand_num <= mut_probability
                mut_cross_mat(r,j) = normrnd(pop_mat(rand_coordinate,j),abs(pop_mat(rand_coordinate,j))*shrink_vals(j));
            else
                mut_cross_mat(r,j) = pop_mat(rand_coordinate,j);
            end
        end
    end
    %% ELITE POPULATION
        % This step generates the elite individuals and places them into the
        % population matrix.
    if nlls_size > 0
        elite_and_nlls_mat(1:elite,1:end) = best_pop(1:elite,1:end);
            % This step introduces the sub-population that was further optimized by
            % NLLS back into the general population.
        elite_and_nlls_mat(elite+1:nlls_size+elite,1:end) = NLLS_pop;
    else
        elite_and_nlls_mat(1:end,1:end) = best_pop(1:elite,1:end);
    end
        % This step combines the elite, nlls, mutation, and crossover populations
    population_mat(1:pop_size-mut_cross_size,1:end) = elite_and_nlls_mat;
    population_mat(pop_size-mut_cross_size+1:pop_size,1:end) = mut_cross_mat;
        %% BOUNDS ON PARAMETERS:
        % The following loop imposes the bounds on each parameter value.
    for k = 1:pop_size
        for i = 1:(size(CyMat,2)-1)
             if population_mat(k,i) >= params_ub(i)
                 population_mat(k,i) = params_ub(i);
             elseif population_mat(k,i) <= params_lb(i)
                 population_mat(k,i) = params_lb(i);
             else
                 population_mat(k,i) = population_mat(k,i);
             end
        end
    end  
    
        %% ITERATION CONTROL:
        % Utilizes iteration control vectors to set which parameters will be
        % fixed or floated in the generations for a given cycle
    for k = 1:pop_size
        for i = 1:(size(CyMat,2)-1)
             if CyMat(i) == 0
                 population_mat(k,i) = params(i);
             else
                 population_mat(k,i) = population_mat(k,i);
             end
        end
    end       
    
    
        %% Insert Previous Best Fit into Population
    if g>1
        load('CurrentParamVals')
        population_mat(end,:)= param_vals;
    end
     %% Recalculation of CHI-SQUARED AND SORTING OF CHI-SQUARED   
        % Chisq values are calculated in parallel in the following parfor loop.
    parfor q = 1:pop_size
        chisq(q,1) = sum(sum(((data-model_simulator(population_mat(q,1:end)',indvar,flag))./weights).^2));     
    end

    [Y,index] = sort(chisq);
    best_chisq = Y(1);
    param_vals = population_mat(index(1),1:end); 
    %% END OF GENERATION DISPLAY:
        % This reports the best Chisq value and the best fit parameter values.
    residual = data-model_simulator(param_vals',indvar,flag);
    RMSD = sqrt(sum(sum(residual.^2))/numel(data));
    VAR = sum(sum(residual.^2))/numel(data);
    simvals = model_simulator(param_vals',indvar1_sim,flag_sim);
        % Creates a matrix of each RMSD,VAR, mean and Chi for every generation
    running_chi(g,1) = Y(1);
    %% PLOTTING AND RUNNING WINDOW
        %This section of code is used to create a user display that shows
        %the experimental data and current simulated data plotted together,
        %a plot of chisq as a function of generation, and a listing of the
        %current and one previous values for all parameters. It also
        %displays which parameters are currently being floated(parametrs
        %turn red), if LSQNL is being used, the current chisq, and the
        %current and total generation, cycle and period number. This was
        %written to specifically for the model used in the construction of the
        %algorithm and will need to be modified in order to fit other
        %models.
    switch run_fig
        case 'on'
            % Creates a growing matrix of generation #
        x_ax = 1:1:g;
        % Figure #1
        runfig = figure(1);
        % Plot Chi-Squared
        subplot(2,2,3)
        semilogy(x_ax,running_chi(1:g,1),'Color','b')
        xlabel('Iteration Number');
        ylabel('Chi-Square');
        title('Chi-Square');
        ytickformat('%.g')
        set(gca,'fontsize',20,'linewidth',2);
        switch axspace
            case 'log'
            % log sapce
        subplot(2,2,[1,2])
        semilogx(indvar,data,'o')
        hold on
        semilogx(indvar1_sim,simvals,'.')
        title('Data/Fits');
        ylabel('Signal')
        set(gca,'fontsize',20,'linewidth',2);
        set(gca,'fontsize',20,'linewidth',2);
        hold off
            case 'lin'
            % lin space
        subplot(2,2,[1,2])
        plot(indvar,data,'o')
        hold on
        plot(indvar1_sim,simvals,'.')
        title('Data/Fits');
        ylabel('Signal')
        set(gca,'fontsize',20,'linewidth',2);
        set(gca,'fontsize',20,'linewidth',2);
        hold off
        end
        % Status and Parameters
        subplot(2,2,4)
        delete(findobj(gca,'Type','Text')) % clears text in each iteration
        o = [153 38 0]./255; % defines the color orange for indicating that params are floating
        text(0,1,'Parameters:','Units','normalized','FontSize',14);
        for i = 1:length(param_vals) % general loop for printing param values
            txt = sprintf('Param %d = %.3f',i,param_vals(i)); % creates text with param inputs
            if CyMat(1) == 1 % if floating disp param in orange
            text(0,1-0.1*i,txt,'Units','normalized','FontSize',14,'Color',o);
            else % if not floating disp param in black
            text(0,1-0.1*i,txt,'Units','normalized','FontSize',14);    
            end
        end
        if CyMat(end) == 1 % if lsqnl is being used disp in orange
            text(0,1.1,'LSQNL ON','Units','normalized','FontSize',14,'Color',o)
        else % if lsqnl is not being used disp in black
            text(0,1.1,'LSQNL OFF','Units','normalized','FontSize',14)
        end
            % disp current chisq
        text(0.5,1.1,['ChiSq: ',sprintf('%.2f',Y(1))],'Units','normalized','FontSize',14,'Color','b')
            % disp current period
        text(0.5,1,['Period: ',s_y],'Units','normalized','FontSize',14)
            % disp current and max cycle
        text(0.5,0.9,['Cycle: ',c_y,' Out of: ',num2str(CyMatH)],'Units','normalized','FontSize',14)
            % disp current and max generation
        text(0.5,0.8,['Generation: ',g_y,' Out of: ',num2str(gzy)],'Units','normalized','FontSize',14)
            % disp current RMSD
        text(0.5,0.7,['RMSD = ',num2str(RMSD)],'Units','normalized','FontSize',14)
            % disp current VAR
        text(0.5,0.6,['VAR = ',num2str(VAR)],'Units','normalized','FontSize',14)
        axis off 
        case 'off'
    end
        drawnow % updates display
%% OUTPUTS     
        %reports the best Chisq in command window and saves the corresponding params
            %and chisq
        chisq_y = num2str(Y(1));    
        disp(['Chisq: ',chisq_y]) %displays best chisq
        disp('  ')
        chi_values(1:end,g)=Y;
    %% PARAMETER SAVE
        %Saves parameters at end of each generation
        save('CurrentParamVals','param_vals');
        saveas(runfig,sprintf('RunFig'));
end 

