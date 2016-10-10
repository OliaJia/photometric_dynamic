%Orthogonal Matching Pursuit to generate c
function [Cs,i,k] = Orthogonal_Matching_Pursuit(A, b, kmax, err)
k = 0; q = b;

% Anor is to store the support vector
Anor=zeros(size(A,1),kmax);

% i is an array that store the index of suport vector
i=zeros(size(A,2),1);


while( norm(q,2)>=err &&k<kmax)
    %find the vector with smallest angle with residual
   [x, i(k+1,1)] = max(abs(q'*(A))); 
   
   %add this vector to the support vector
    Anor(:,k+1)=A(:,i(k+1,1));
    
    % find the coefficient of this support vector with least square
    cs = Anor(:,1:k+1)\b;
    
    %calculate the residual
    q = b - Anor(:,1:k+1)*cs;
    
    k=k+1;
end
Cs=zeros(size(A,2),1);

% place the coefficient to the right place in the similarity matrix.
for l=1:1:k
    Cs(i(l,1),1)=cs(l,1);
end
    


  
