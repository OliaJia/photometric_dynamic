function group_real=recover_real_group(N,w,groups)
group_real=zeros(N,1);
z=zeros(size(w,1),1);
w_new=[w z];
for i=1:1:size(groups,1)
    m=1;
    while(w_new(i,m)~=0)
        group_real(w_new(i,m),1)=groups(i,1);
        m=m+1;
    end
end
