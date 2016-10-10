function [X,Xs,w]=Matrix_Collapse_2(X,ep)
n=size(X,2);
ind=1;
%% collapse the matrix with greedy search method
% % 
T=6400;
for m=1:1:12
    for i=1:1:T
        if  sum(X(:,(m-1)*T+i),1) ~= 0
            w(ind,1) =(m-1)*T+i;
            p = 2;
            for j = i+1 : 1 : T
                if sum(X(:,(m-1)*T+j),1) == 0
                    continue;
                else
                    t=sum(abs(X(:,(m-1)*T+i)-X(:,(m-1)*T+j)),1);
                    if  t< ep
                        X(:,(m-1)*T+j)=0;
                        w(ind,p)=(m-1)*T+j;
                        p=p+1;
                    end
                end
            end
            ind=ind + 1;
        end
    end
end

%% generate the collapse matrix
m=1;Xs=zeros(6,3);
for i=1:1:n
    if (sum(X(:,i),1)~=0)
        Xs(:,m)=X(:,i);
        m=m+1;
    end
end
