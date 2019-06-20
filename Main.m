clc
clear 
close all

%% set paths

addpath(genpath('functions'))
my_dir = pwd; 
% parpool(16)

%% variables initialization 

noise_level = [0,0.5,1];% 0:0.1:1; 
num_test = 3; 
num_trials = 30; 

acc = cell(1,num_test); 
TE_cell = cell(1,num_test);
myswitch = zeros(1,num_test); 

for k = 1 : num_test

    myswitch(k) = randn(1); % Directionality parameter
    acc_aux = zeros(numel(noise_level),num_trials); 
    TE_aux = cell(numel(noise_level),1); 

    iter = 1;

    for noise = noise_level            
        %% data generation 
        fprintf('\n Test %d/%d, Noise level: %.2f \n',k,num_test,noise)

        % sampling frequency 
        fsample = 256; 

        % labels 
        label = {'X','Y'};

        % trials data 
        N = 512;    % Time points 
        P = 3;      % Model order
        trial = cell(1,num_trials); 
        time = cell(1,num_trials); 

        for i=1:num_trials
            stability_check = 0;
            while stability_check == 0
                trial_aux = simuldata_uniform(N,P,noise,myswitch(k));
                stability_check = (abs(trial_aux(end,1)/trial_aux(1,1))<3)*(abs(trial_aux(end,2)/trial_aux(1,2))<3);
            end 
            trial{i} = trial_aux'; 
            time{i} = (0:size(trial_aux,1)-1)/fsample;
        end

        % data
        data.fsample=fsample;
        data.label=label;
        data.trial=trial;
        data.time=time; 

        %% Calculating TE

        fprintf('\n Computing TE ... \n')
        dim = 3;    %embedding dimension
        tau_aux = 1; 
        u = 1; 
        maxlag = 20; 
        
        TE_1 = zeros(num_trials,1);
        TE_2 = zeros(num_trials,1);

        for i = 1:num_trials
            fprintf('\n Trial %d/%d',i,num_trials)

            ts_1 = data.trial{1,i}(1,:);
            ts_2 = data.trial{1,i}(2,:);

            alpha=2;  

            % X -> Y 
            ACT = ACT_estimation(ts_2,maxlag);
            tau = round(ACT*tau_aux);  %embedding delay in number of sampled points
            teS_1 = kRTE(ts_1,ts_2,dim,tau,u,alpha,1);
            TE_1(i) = teS_1; 

            % Y -> X 
            ACT = ACT_estimation(ts_1,maxlag);
            tau = round(ACT*tau_aux);  %embedding delay in number of sampled points
            teS_2 = kRTE(ts_2,ts_1,dim,tau,u,alpha,1);
            TE_2(i) = teS_2;

        end    
        TE = [TE_1,TE_2];
        
    %% Accuracy 

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

% Accuracy plots 

acc_aux = acc(1,:); 
acc_var = reshape(mean(cat(3,acc_aux{:}),2),numel(noise_level),num_test);

figure
lineProps ={'-k','Linewidth',2};
H1=shadedErrorBar(noise_level,acc_var', {@mean, @(features) 1*std(features)  }, lineProps, 0.5);
ylim([0,100]), grid on 
xlabel('Noise level')
ylabel('Accuracy (%)')

legend([H1.patch],'kRTE \alpha=2','location','southeast')



