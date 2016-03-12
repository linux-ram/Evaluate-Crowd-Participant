function [p_ann,score1]=adversarial(Pz_x_b,Py_xzbar,Py_xz,number_of_labelers)

p_ann=[]; p_ann_temp1=0; p_ann_temp2=0;

for i=1:size(Py_xzbar,1)
    for j=1:size(Py_xzbar,2)
        
        p_ann_temp1=(1-Pz_x_b(i))*prod(Py_xzbar(i,:))+Pz_x_b(i)*prod(Py_xz(i,:)); %numerator in eq. 19
        
        Pyxzbar=Py_xzbar(i,:); Pyxz=Py_xz(i,:);
        
        Pyxzbar(j)=1; Pyxz(j)=1;
        
        p_ann_temp2=(1-Pz_x_b(i))*prod(Pyxzbar)+Pz_x_b(i)*prod(Pyxz); %denominator in eq. 19
        p_ann(i,j)=p_ann_temp1/p_ann_temp2; %The term on the left in Eq. 19 computed for every annotator and datapoint
        
    end
end

%computing scores for each annotator
%product of p_ann's corresponding to each annotator
score1=zeros(1,number_of_labelers);
for j=1:number_of_labelers
    p_ann(isnan(p_ann(:,j)),j)=1;
    score1(j)=-sum(log(p_ann(:,j))); %negative sum of log's of p_ann for every annotator
end

%
% figure
% plot([1:size(p_ann,1)],p_ann(:,1),'+',[1:size(p_ann,1)],p_ann(:,2),'*',[1:size(p_ann,1)],p_ann(:,3),'o')
% axis([1,size(p_ann,1),0,1])
% title('Evaluating Annotators: Eq. 19 for each datapoint [no adversaries, binomial case, random initialization]')
% xlabel('Datapoint');
% ylabel('y^{k}/y^{(t\\k)},x');
% %legend('Annotator 1','Annotator 2','Annotator 3');
% legend(['Annotator 1 Score: ' num2str(score(1))],['Annotator 2 Score: ' num2str(score(2))], ['Annotator 3 Score: ' num2str(score(3))]);


% %breast dataset
% %computing scores for each annotator
% score=[0 0 0 0 0];
% for j=1:5
% score(j)=-sum(log(p_ann(:,j)));
% end
%
% figure
% plot([1:size(p_ann,1)],p_ann(:,1),'+',[1:size(p_ann,1)],p_ann(:,2),'*',[1:size(p_ann,1)],p_ann(:,3),'o',[1:size(p_ann,1)],p_ann(:,4),'s',[1:size(p_ann,1)],p_ann(:,5),'d')
% axis([1,size(p_ann,1),0,1])
% title('Evaluating Annotators: Eq. 19 for each datapoint [Breast Dataset, 3 annotators and two adversaries with 0.2 flip prob. resp., binomial case, random initialization]')
% %title('Evaluating Annotators: Eq. 19 for each datapoint [ Dataset,adversary: Annotator 3 with p=0.6, binomial case, random initialization]')
% xlabel('Datapoint');
% ylabel('y^{k}/y^{(t\\k)},x');
% %legend('Annotator 1','Annotator 2','Annotator 3');
%
% % k=5; a=[];
% % for j=1:k
% %     if j~=k
% %         a=[[a,['[Annotator ' num2str(j) ' Score: ' num2str(score(j))]] '], '];
% %     else
% %         a=[a,['[Annotator ' num2str(j) ' Score: ' num2str(score(j))] ']'];
% %     end
% % end
% legend(['Annotator 1 Score: ' num2str(score(1))],['Annotator 2 Score: ' num2str(score(2))],['Annotator 3 Score: ' num2str(score(3))],['Adversary 1 Score: ' num2str(score(4))],['Adversary 2 Score: ' num2str(score(5))]);
% % legend(a);
