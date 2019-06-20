function [TE,K,H] = kRTE(X,Y,dim,tau,u,alpha,var)
    
    % Kernel transfer entropy 
    
    T = length(X);          % length of full time series
    L = T-((dim-1)*tau);    % number of points inside the time series available for advance and delay embedding
    FirstP = T-L;           % all referencing of embeddings is done wrt the prediction points

    Y_emb = zeros(L,dim);
    X_emb = zeros(L,dim); 

    for ii = 1:L           
        for jj = 1:dim
            Y_emb(ii,jj) = Y(ii+FirstP-(jj-1)*tau); 
            X_emb(ii,jj) = X(ii+FirstP-(jj-1)*tau);
        end
    end

    Y_t = Y(FirstP+u+1:end)'; 
    Y_emb = Y_emb(u:end-1,:);
    X_emb = X_emb(1:end-u,:);

    sigma = var*median(pdist(Y_emb)); %kScaleOptimization_info(pdist(Y_emb_t));
    K_Y_emb = exp(-(pdist2(Y_emb,Y_emb).^2)/(2*sigma^2));

    sigma = var*median(pdist(X_emb));
    K_X_emb = exp(-(pdist2(X_emb,X_emb).^2)/(2*sigma^2));

    sigma = var*median(pdist(Y_t));
    K_Y_t = exp(-(pdist2(Y_t,Y_t).^2)/(2*sigma^2));

    K1 = (K_Y_emb.*K_X_emb)/trace(K_Y_emb.*K_X_emb); 
    H1 = kernel_entropy(K1,alpha); 

    K2 = (K_Y_t.*K_Y_emb.*K_X_emb)/trace(K_Y_t.*K_Y_emb.*K_X_emb);
    H2 = kernel_entropy(K2,alpha); 

    K3 = (K_Y_t.*K_Y_emb)/trace(K_Y_t.*K_Y_emb); 
    H3 = kernel_entropy(K3,alpha); 

    K4 = K_Y_emb/trace(K_Y_emb); 
    H4 = kernel_entropy(K4,alpha); 

    K = {K1,K2,K3,K4};
    H = {H1,H2,H3,H4};
    
    TE = H1-H2+H3-H4; 
end 