function H = kernel_entropy(A,alpha)
    H = real((1/(1-alpha))*log2(trace(A^alpha))); 
end 