function H = kernel_entropy(A,alpha)

    % Kernel-based Renyi's entropy
    
    % Inputs:
    % A in R^{N x N} : Input kernel matrix 
    % alpha : Renyi's entropy order 
    
    % Output:
    % H : Kernel-based Renyi entropy
    
    % Ivan De La Pava Panche, Automatics Research Group
    % Universidad Tecnologica de Pereira, Pereira - Colombia
    % email: ide@utp.edu.co

    H = real((1/(1-alpha))*log2(trace(A^alpha))); 
end 