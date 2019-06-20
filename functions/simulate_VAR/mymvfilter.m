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