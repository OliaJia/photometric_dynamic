m=1;Xs=zeros(6,3);
for i=1:1:n
    if (X(:,i)~=0)
        Xs(:,m)=X(:,i);
        m=m+1;
    end
end
