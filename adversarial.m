function [p_ann,score1]=adversarial(Pz_x_b,Py_xzbar,Py_xz,number_of_labelers)

p_ann=[]; p_ann_temp1=0; p_ann_temp2=0;

for i=1:size(Py_xzbar,1)
    for j=1:size(Py_xzbar,2)
        %numerator in eq. 19
        p_ann_temp1=(1-Pz_x_b(i))*prod(Py_xzbar(i,:))+Pz_x_b(i)*prod(Py_xz(i,:));
        Pyxzbar=Py_xzbar(i,:); Pyxz=Py_xz(i,:);
        Pyxzbar(j)=1; Pyxz(j)=1;
        %denominator in eq. 19
        p_ann_temp2=(1-Pz_x_b(i))*prod(Pyxzbar)+Pz_x_b(i)*prod(Pyxz);
        %The term on the left in Eq. 19 computed for every annotator and datapoint
        p_ann(i,j)=p_ann_temp1/p_ann_temp2;
    end
end
%computing scores for each annotator
%product of p_ann's corresponding to each annotator
score1=zeros(1,number_of_labelers);
for j=1:number_of_labelers
    p_ann(isnan(p_ann(:,j)),j)=1;
    %negative sum of log's of p_ann for every annotator
    score1(j)=-sum(log(p_ann(:,j)));
end