X=Xnor;ep=0.00001;
n=size(X,2);
ind=1;
for i=1:1:n
    if  sum(X(:,i),1) ~= 0
        w(ind,1) = i; 
        m = 2;
        for j = i+1 : 1 : n
            if sum(X(:,j),1) == 0
                continue;
            else
                if norm(X(:,i)-X(:,j),2) < ep
                    X(:,j)=0;
                    w(ind,m)=j;
                    m=m+1;
                end
            end
        end
        ind=ind + 1;
    end
end
m=1;Xs=zeros(6,3);
for i=1:1:n
    if (sum(X(:,i),1)~=0)
        Xs(:,m)=X(:,i);
        m=m+1;
    end
end
