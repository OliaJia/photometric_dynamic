function  [Cs,i,k] =Nearest_Subspace_Neighbor(A, b, kmax, d, num)
i=zeros(size(A,2),1);
Cs=zeros(size(A,2),1);Cs(num,1)=0;
Anor=zeros(size(A,1),1); Anor(:,1)=b;
for k=1:1:kmax
    if  k<=d
        [u,~,~]=svd(Anor,0);
    end
    temp=u'*A;
    [~, i(k,1)] = max(sum(temp.*temp,1)); 
    Anor(:,k+1)=A(:,i(k,1)); A(:,i(k,1))=0;
   % Cs(i(k),1)=0.95^(k);
end
   cs = Anor(:,2:k+1)\b; 
for l=2:1:k
    Cs(i(l,1),1)=cs(l,1);
end
