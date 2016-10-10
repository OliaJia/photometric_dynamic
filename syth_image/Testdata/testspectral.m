C_spectral_1=zeros(size(C,2),size(C,2));
for i=1:1:size(C,2)
    for j=1:1:size(C,1)
        C_spectral_1(i,floor(Id(j,i)))=abs(C(j,i));
    end
end
C_spectral=C_spectral_1'+C_spectral_1;
groups = SpectralClustering(C_spectral,3);
