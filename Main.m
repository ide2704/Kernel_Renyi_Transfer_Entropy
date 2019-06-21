% Kernel-based Renyi transfer entropy

% This script generates the results for the kernel-based Renyi transfer 
% entropy method (alpha = 2) shown in Figure 2A of the paper "A data-driven 
% measure of effective connectivity based on Renyi’s alpha-entropy". 

% Ivan De La Pava Panche, Automatics Research Group
% Universidad Tecnologica de Pereira, Pereira - Colombia
% email: ide@utp.edu.co

clc
clear 
close all
    
% Set paths

addpath(genpath('functions'))
my_dir = pwd; 

% Variables initialization 

noise_level = 0:0.1:1; % Noise levels
num_test = 10;         % Number of realizations 
num_trials = 100;      % Numebr of trials per realization (For the script to finish faster less trials can be selected) 

acc = cell(1,num_test); 
TE_cell = cell(1,num_test);
myswitch = zeros(1,num_test); 

for k = 1 : num_test

    myswitch(k) = randn(1); % Directionality parameter
    acc_aux = zeros(numel(noise_level),num_trials); 
    TE_aux = cell(numel(noise_level),1); 

    iter = 1;

    for noise = noise_level            
        % Data generation 
        fprintf('\n Test %d/%d, Noise level: %.2f \n',k,num_test,noise)

        % Trials data 
        N = 512;    % Time points (data size)
        P = 3;      % Autoregressive model order
        trial = cell(1,num_trials); 

        for i=1:num_trials
            stability_check = 0;
            while stability_check == 0
                trial_aux = simuldata_uniform(N,P,noise,myswitch(k));
                stability_check = (abs(trial_aux(end,1)/trial_aux(1,1))<3)*(abs(trial_aux(end,2)/trial_aux(1,2))<3);
            end 
            trial{i} = trial_aux'; 
        end
        data.trial=trial;

        % Kernel-based Renyi transfer entropy estimation 
        fprintf('\n Computing TE ... \n')
        
        dim = 3;     % Embedding dimension 
        tau_aux = 1; % Embedding delay (in ACT)
        u = 1;       % Interaction time (in samples)
        maxlag = 20; 
    
        TE_1 = zeros(num_trials,1);
        TE_2 = zeros(num_trials,1);

        for i = 1:num_trials
            fprintf('\n Trial %d/%d',i,num_trials)

            ts_1 = data.trial{1,i}(1,:);
            ts_2 = data.trial{1,i}(2,:);

            alpha = 2; % Renyi's entropy order  

            % TE (X -> Y) 
            ACT = ACT_estimation(ts_2,maxlag); % Autocorrelation time 
            tau = round(ACT*tau_aux);          % Embedding delay in number of sampled points
            TE_1(i) = kRTE(ts_1,ts_2,dim,tau,u,alpha); % Transfer entropy estimation  

            % TE (Y -> X) 
            ACT = ACT_estimation(ts_1,maxlag); % Autocorrelation time
            tau = round(ACT*tau_aux);          % Embedding delay in number of sampled points
            TE_2(i) = kRTE(ts_2,ts_1,dim,tau,u,alpha); % Transfer entropy estimation

        end    
        TE = [TE_1,TE_2]; % TE matrix 
        
    % Accuracy 
        if myswitch(k) > 0
            acc_aux(iter,:) = 100*((TE(:,1)-TE(:,2))>0);
        else
            acc_aux(iter,:) = 100*((TE(:,2)-TE(:,1))>0);
        end

        fprintf('\n\n Mean TE accuracy: %.2f, iter %d\n',mean(acc_aux(iter,:)),iter)

        TE_aux{iter} = TE; 
        iter = iter+1; 
    end 
    
    acc{k} = acc_aux; 
    TE_cell{k} = TE_aux;
end 

% Accuracy plot

acc_aux = acc(1,:); 
acc_var = reshape(mean(cat(3,acc_aux{:}),2),numel(noise_level),num_test);

figure
lineProps ={'-k','Linewidth',2};
H1=shadedErrorBar(noise_level,acc_var', {@mean, @(features) 1*std(features)  }, lineProps, 0.5);
ylim([0,100]), grid on 
xlabel('Noise level')
ylabel('Accuracy (%)')

legend([H1.patch],'kRTE \alpha=2','location','southeast')



