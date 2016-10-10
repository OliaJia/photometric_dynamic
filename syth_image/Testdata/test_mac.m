for i=1:1:size(w,1)
    for j=1:1:size(w,2)
        if(w(i,j)==0)
            continue;
        else
            w(i,j)=img(floor(w(i,j)),1);
        end
    end
end
