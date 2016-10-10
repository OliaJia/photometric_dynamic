function Cr=Matrix_recover(C,Id)

N=size(C,2);
Cr=zeros(N,N);
for i=1:1:N
    for j=1:1:size(C,1)
        Cr(i,Id(j,i))=C(j,i);
    end
end

