function Omega = Contr_matrix(nX,Sel,class)

Omega=zeros(nX,nX);
for k=1:class-1
    a=Sel{k};b=[];
    for r=k+1:class
        b=[b;Sel{r}(:)];
    end
    
    for i=1:length(a)
        for j=1:length(b)
            Omega(a(i),b(j))=1;
        end
    end
end
Omega=Omega +Omega' +eye(nX); 