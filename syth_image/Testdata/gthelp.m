num=zeros(26,2);
for i=1:1:26
    num(i,1)=size(find(img==i),1);
    num(i,2)=size(find(group_real==i),1);
end
Repre1=zeros(26,7);
for i=1:1:26
    a=find(img==i);
    Repre1(i,:)=code_point(a(1,1),:);
    %clear a
end
Repre2=zeros(26,7);
for i=1:1:26
    a=find(group_real==i);
    Repre2(i,:)=code_point(a(1,1),:);
    %clear a
end
