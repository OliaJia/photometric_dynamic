%% Collapse the matrix, if the the two vector are really close to each other, we just treat them as one
% % we do this blockly , which means this method may may not collapse all
% the vector near to each other
function [X,Xs,w]=Matrix_Collapse(X,ep)
n=size(X,2);
ind=1;
%% collapse the matrix with greedy search method
% % 
for i=1:1:n
        if  sum(X(:,i),1) ~= 0
            w(ind,1) =i;
            p = 2;
            for j = i+1 : 1 : n
                if sum(X(:,j),1) == 0
                    continue;
                else
                    t=sum(abs(X(:,i)-X(:,j)),1);
                    if  t< ep
                        X(:,j)=0;
                        w(ind,p)=j;
                        p=p+1;
                    end
                end
            end
            ind=ind + 1;
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
