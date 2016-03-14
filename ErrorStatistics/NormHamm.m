function [acc det_rat fa_rat]=NormHamm(tg, reg)

%%...this is the function that is for checking out the accuracy 
%%...between the real target and regression
%%...also is known as Hamming distance

%%...codes starts here
acc=0;
det_rat=0;
fa_rat=0;
pos_num=0;
neg_num=0;
if size(tg, 1)<size(tg, 2)
    tg=tg';
    reg=reg';
end
for i=1:1:size(tg,1)
    %if tg(i,1)==reg(i,1)  %%ver 1.
    if tg(i,1)*reg(i,1)>0
        acc=acc+1;
    end
    
    if tg(i,1)>0
        pos_num=pos_num+1;
        if reg(i,1)>0
            det_rat=det_rat+1;
        end
    end
    
    if tg(i,1)<=0
        neg_num=neg_num+1;
        if reg(i,1)>0
            fa_rat=fa_rat+1;
        end
    end
end
acc=acc/size(tg, 1);
if pos_num>0
    det_rat=det_rat/pos_num;
else
    det_rat=1;
end
if neg_num>0
    fa_rat=fa_rat/neg_num;
else
    fa_rat=0;
end

return