function TE = kRTE(X,Y,dim,tau,u,alpha)
    
    % Kernel-based Renyi transfer entropy
    
    % Inputs:
    % X : Driving time series (in R^{N})
    % Y : Driven time series (in R^{N})
    % dim : Embedding dimension 
    % tau : Embedding delay (in samples)
    % u : Interection time (in samples)
    % alpha : Renyi's entropy order 
    
    % Output:
    % TE :  Kernel-based Renyi transfer entropy
    
    % Ivan De La Pava Panche, Automatics Research Group
    % Universidad Tecnologica de Pereira, Pereira - Colombia
    % email: ide@utp.edu.co
        
    % Time embedding 
    T = length(X);          % Length of full time series
    L = T-((dim-1)*tau);    % Number of points inside the time series available for advance and delay embedding
    FirstP = T-L;           % All referencing of embeddings is done with regard to the prediction points
    
    Y_emb = zeros(L,dim);
    X_emb = zeros(L,dim); 

    for ii = 1:L           
        for jj = 1:dim
            Y_emb(ii,jj) = Y(ii+FirstP-(jj-1)*tau); 
            X_emb(ii,jj) = X(ii+FirstP-(jj-1)*tau);
        end
    end
    
    % Temporal adjustment of the embedded time series according to the interaction time u 
    Y_t = Y(FirstP+u+1:end)'; 
    Y_emb = Y_emb(u:end-1,:);
    X_emb = X_emb(1:end-u,:);

    % Kernel matrices computation from the embedded time series 
    sigma = median(pdist(Y_emb)); % Kernel bandwidth
    K_Y_emb = exp(-(pdist2(Y_emb,Y_emb).^2)/(2*sigma^2));

    sigma = median(pdist(X_emb));
    K_X_emb = exp(-(pdist2(X_emb,X_emb).^2)/(2*sigma^2));

    sigma = median(pdist(Y_t));
    K_Y_t = exp(-(pdist2(Y_t,Y_t).^2)/(2*sigma^2));

    % Joint and marginal kernel-based Renyi entropy estimation  
    K1 = (K_Y_emb.*K_X_emb)/trace(K_Y_emb.*K_X_emb); 
    H1 = kernel_entropy(K1,alpha); 

    K2 = (K_Y_t.*K_Y_emb.*K_X_emb)/trace(K_Y_t.*K_Y_emb.*K_X_emb);
    H2 = kernel_entropy(K2,alpha); 

    K3 = (K_Y_t.*K_Y_emb)/trace(K_Y_t.*K_Y_emb); 
    H3 = kernel_entropy(K3,alpha); 

    K4 = K_Y_emb/trace(K_Y_emb); 
    H4 = kernel_entropy(K4,alpha); 
    
    % Transfer entropy estimation (sum of kenel-based Renyi entropies) 
    TE = H1-H2+H3-H4; 
end 