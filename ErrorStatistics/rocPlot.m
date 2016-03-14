function [faAxe detAxe]=rocPlot(TrueLab, Regression, dispOpt, Smooth, evenAxe)
if nargin<5
    evenAxe=0;
end
if nargin<4
    Smooth=0;
end
if nargin<3
    dispOpt=0;
end
if Smooth==1
    Regression=Regression+0.01*randn(size(Regression,1),1);
end
Bin=sort(unique(Regression), 'descend')';
faAxe=[];
detAxe=[];
for i=1:1:size(Bin,2)
    CurRegCant=zeros(size(Regression,1),1);
    for j=1:1:size(Regression, 1)
        if Regression(j,1)>=Bin(1,i)
            CurRegCant(j,1)=1;
        end
    end
    
    [acc curdet curfa]=NormHamm(TrueLab, CurRegCant);
    faAxe=[faAxe curfa];
    detAxe=[detAxe curdet];
end
faAxe=[0 faAxe 1];
detAxe=[0 detAxe 1];


if evenAxe==1
    xAxe=[0:0.025:1];    
    yAxe=zeros(1,size(xAxe,2));
    for i=1:1:size(xAxe,2)
        startSet=find(faAxe<=xAxe(1,i));
        startPoint=startSet(1,end);
        endSet=find(faAxe>=xAxe(1,i));
        endPoint=startSet(1,end);
        yAxe(1,i)=(detAxe(1,startPoint)+detAxe(1,endPoint))/2;
    end
    faAxe=xAxe;
    detAxe=yAxe;
end

if dispOpt==1
    plot(faAxe, detAxe, '--'), grid on;
end

return