function [Pz Z]=MjVt(Y)

[n t]=size(Y);
Pz=zeros(n,1);
Z=zeros(n,1);
for i=1:1:n
    ctP=0;
    ctN=0;
    for j=1:1:t
        if Y(i,j)==1
            ctP=ctP+1;
        elseif Y(i,j)==-1
            ctN=ctN+1;
        end
    end
    if ctP==ctN
        Pz(i,1)=0.5;
        Z(i,1)=1;
    elseif ctP==0
        Pz(i,1)=0.1;
        Z(i,1)=-1;
    elseif ctN==0
        Pz(i,1)=0.9;
        Z(i,1)=1;
    else
        Pz(i,1)=min(0.9, ctP/(ctP+ctN));
        if Pz(i,1)>=0.5
            Z(i,1)=1;
        else
            Z(i,1)=-1;
        end
    end
end
        
return