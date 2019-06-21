function [data]=simuldata_uniform(N,P,gamma,myswitch)
% Causality Challenge (Estimate causal direction in noisy background)
% http://clopinet.com/causality/data/nolte/

%generates one example for the causality challenge
%
% Input:
% N : number of time points
% P : order of AR-system
% gamma: parameter controlling the relative strength of noise 
%        and signal, gamma=0 means only signal, gamma=1 means only noise
% myswitch:  if positive direction is from channel 1 to channel 2;
%            if negative direction is from channel 2 to channel 1;
%
% Output: 
% data   Nx2 matrix of data
%
% One example of the challenge was generated with the command:
% N=6000; P=10; gamma=rand;myswitch=randn;
% data=simuldata(N,P,gamma,myswitch)
% 


ssig=1-gamma;
snoise=gamma;
nrun=1;
M=2;M1=2;
M2=3;
ckount=0;
while ckount< nrun
    Arsig=[];
    for k=1:P
      aloc=randn(M)/4;
      if myswitch<0
          aloc(2,1)=0;
%           aloc(1,2)=1;
      else
          aloc(1,2)=0;
%           aloc(2,1)=1;
      end
      Arsig=[Arsig,aloc];
    end
    E=eye(M*P);AA=[Arsig;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));

    Arnoise=[];
    for k=1:P;aloc=diag(diag(randn(M2)/4)); Arnoise=[Arnoise,aloc];end
    
    E=eye(M2*P);AAnoise=[Arnoise;E(1:end-M2,:)];
    lambda=eig(AAnoise);
    lambdamaxnoise=max(abs(lambda));

    %subplot(3,3,i);plot(imag(ps));
    if lambdamax<1 && lambdamaxnoise<1
        
        ckount=ckount+1;
         x=rand(M1,N)-.5; 
         %x=randn(M1,N); 
         data= mymvfilter(Arsig,x);
         datasig=data';
         siglevel=norm(datasig,'fro');
         x=rand(M2,N)-.5; 
         %x=randn(M2,N); 
         data= mymvfilter(Arnoise,x);
         datanoise=data';
         B=randn(M2,M1);datanoise=datanoise*B;
         noiselevel=norm(datanoise,'fro');
         data=ssig*datasig/siglevel+snoise*datanoise/noiselevel;
    end
end

return;

function y=mymvfilter(Ar,x)
 
 [nchan,norder]=size(Ar);
 norder=norder/nchan;
 Ar_reshape=zeros(nchan,nchan,norder);
 for i=1:norder
  Ar_reshape(:,:,i)=Ar(:,(i-1)*nchan+1:i*nchan);
 end
 
 [~,N]=size(x);
 y=x;
 for i=2:N
     norderloc=min([i-1,norder]);
     for j=1:norderloc
         y(:,i)=y(:,i)+Ar_reshape(:,:,j)*y(:,i-j);
     end
 end
         
return

% http://clopinet.com/causality/data/nolte/

%To read the data into MATLAB, type 
%fid=fopen('simuldata.bin');
%data=reshape(fread(fid,'float'),6000,2,1000); 

%To read the data e.g. of the first subject into Matlab type:
%fid=fopen('sub1.bin','r');
%data=reshape(fread(fid,'float'),[],19)




