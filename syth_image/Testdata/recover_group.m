%recover the real cluster from the C
group_real=zeros(size(Xnor,2),1);
z=zeros(size(w,1),1);
w_new=[w z];
for i=1:1:size(groups,1)
    m=1;
    while(w_new(i,m)~=0)
        group_real(w_new(i,m),1)=groups(i,1);
        m=m+1;
    end
end
