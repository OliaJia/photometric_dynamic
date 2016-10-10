function [best_out,best_in]=light_increment(X_out,lt_sel,lt_left,NL)
[U,S,V]=svd(X_out,0);

US=U(:,1:NL)*S(1:NL,1:NL)';
best_out=zeros(size(lt_left,2),1);
best_in=zeros(size(lt_left,2),2);
ind_sample=zeros(3,2);
for i=1:1:size(lt_left,2)
    [ind_sample(1,:),B_1]=greedy_omp_matrix(X_out(:,lt_sel),US,2);
    [~,B_1]=greedy_omp_matrix(X_out(:,lt_left(i)),B_1,1);
    [ind_sample(2,1),B_2]=greedy_omp_matrix(X_out(:,lt_sel),US,1);
    [~,B_2]=greedy_omp_matrix(X_out(:,lt_left(i)),B_2,1);
    [ind_sample(2,2),B_2]=greedy_omp_matrix(X_out(:,lt_sel),B_2,1);
    [~,B_3]=greedy_omp_matrix(X_out(:,lt_left(i)),US,1);
    [ind_sample(3,:),B_3]=greedy_omp_matrix(X_out(:,lt_sel),B_3,2);
    [best_out(i,1),ind]=min([norm(B_1),norm(B_2),norm(B_3)]);
    best_in(i,:)=ind_sample(ind,:);
end

