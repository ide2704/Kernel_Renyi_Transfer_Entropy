function ACT = ACT_estimation(data,maxlag)
    % Autocorrelation estimation 
    
    % Modified version of the functions TEactdetect and TEautocorr from the 
    % TE toolbox TRENTOOL (Version 2.0 by Raul Vicente, Michael Wibral, 
    % and Michael Lindner, Bonn 2011)
    

    %% Remember the working directory

    thresh = exp(-1);

    data_cut = data; 
    data_corr = data_cut-mean(data_cut);

    %% calculating the autocorrelation of the data
    data = data_corr'; 

    % size of data
    [M,N] = size(data);

    % find nearest power
    m = 2*M-1;
    [F,E] = log2(abs(m));

    % Check if m is an exact power of 2.
    if ~isempty(F) && F == 0.5
        E = E-1;
    end

    % check infinity
    checkinfinity = ~isfinite(F);
    E(checkinfinity) = F(checkinfinity);

    % Compute autocorrelation via FFT
    fft_data = fft(data,2^E);
    c = ifft(abs(fft_data).^2);

    % Force real data
    c = real(c);

    % define lags for moving
    if maxlag >= M
        c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
    else
        c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
    end
    % normalize data
    c = c./c(maxlag+1);

    d = c(maxlag+1:end);
    auxlag = 0:maxlag;

    if ~isempty(min(auxlag(d<thresh)))
        ACT = min(auxlag(d<thresh)); % ACT in samples (!)
    elseif isempty(min(auxlag(d<thresh)))
%         ACT = inf; %undetectable ACT is set to infitiy -> ACT threshhold in function transferentropy ACT<threshold
        ACT = maxlag;
    end
end 