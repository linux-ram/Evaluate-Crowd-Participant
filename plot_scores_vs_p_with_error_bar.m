function plot_scores_vs_p_with_error_bar(mean_score,err_bar,adversaries,p)
%mean score vs flip prob. p with std. error bars for each annotator
%adversaries is an array of indices corresponding to the adversaries. Their annotations usually correspond to the rightmost appended columns.
%flag identifies the type of score computed



number_of_annotators=adversaries(end); % need to modify to run the whole code with no adversaries, i.e. k=0
anno_marker={'k+:' 'bd:' 'go:' 'c*:' 'mv:' 'y+:'};
adv_marker={'r>-' 'rp-' 'r^-' 'rs-' 'rh-' 'rd-' 'r+-' 'ro-' 'rv-' 'r<-'};
% legend_ann={'Annotator 1' 'Annotator 2' 'Annotator 3' 'Annotator 4' 'Annotator 5' 'Annotator 6' 'Annotator 7'};
% legend_adv={'Adversary 1' 'Adversary 2' 'Adversary 3' 'Adversary 4' 'Adversary 5' 'Adversary 6' 'Adversary 7' 'Adversary 8' 'Adversary 9' 'Adversary 10'};
%for loop
%plot(p,mean_score(:,1),'k+:',p,mean_score(:,2),'bd:',p,mean_score(:,3),'go:',p,mean_score(:,4),'r>',p,mean_score(:,5),'rp',p,mean_score(:,6),'rs',p,mean_score(:,7),'rh',p,mean_score(:,8),'rd',p,mean_score(:,9),'r+',p,mean_score(:,10),'r^',p,mean_score(:,11),'ro',p,mean_score(:,12),'r+')
%plot(p,mean_score(:,1),'k+:',p,mean_score(:,2),'bd:',p,mean_score(:,3),'go:',p,mean_score(:,4),'r>',p,mean_score(:,5),'rp',p,mean_score(:,6),'r^')
% plot(p,mean_score(:,1),'k+:',p,mean_score(:,2),'bd:',p,mean_score(:,3),'go:',p,mean_score(:,4),'c*:',p,mean_score(:,5),'mv:',p,mean_score(:,6),'r>',p,mean_score(:,7),'rs',p,mean_score(:,8),'rh',p,mean_score(:,9),'rd',p,mean_score(:,10),'r^',p,mean_score(:,11),'rv',p,mean_score(:,12),'r+',p,mean_score(:,13),'rd',p,mean_score(:,14),'ro',p,mean_score(:,15),'r<')
% plot(p,mean_score(:,1),'k+:',p,mean_score(:,2),'bd:',p,mean_score(:,3),'go:',p,mean_score(:,4),'mv:',p,mean_score(:,5),'r>')

% Plotting score curves for annotators and adversaries
for i=1:adversaries(1)-1
plot(p,mean_score(:,i),anno_marker{i})
hold on
end
clear i
j=1;
for i=adversaries
plot(p,mean_score(:,i),adv_marker{j})
j=j+1;
hold on
end
clear i j
errorbar(repmat(p',1,number_of_annotators),mean_score,err_bar,'kx')
title('Annotator score variation for various p_a')
xlabel('Adversary p_a - flip probability')
ylabel('Score for every annotator')


% for red=adversaries
% plot(p,mean_score(:,red),'r') %plot the adversaries discrete scores interpolated in solid red line
% end

%legend('Annotator 1','Annotator 2','Annotator 3','Adversary 1','Adversary 2','Adversary 3','Adversary 4','Adversary 5','Adversary 6','Adversary 7','Adversary 8','Adversary 9');
%legend('Annotator 1','Annotator 2','Annotator 3','Adversary 1','Adversary 2','Adversary 3','Adversary 4');
%legend('Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1','Adversary 2','Adversary 3','Adversary 4','Adversary 5','Adversary 6','Adversary 7','Adversary 8','Adversary 9','Adversary 10');
%legend('Annotator 1','Annotator 2','Annotator 3','Annotator 4','Adversary 1');

% Legend definitions
legend_ann{3,1}={'Annotator 1','Annotator 2','Annotator 3','Adversary 1'};
legend_ann{3,3}={'Annotator 1','Annotator 2','Annotator 3','Adversary 1','Adversary 2','Adversary 3'};
legend_ann{3,9}={'Annotator 1','Annotator 2','Annotator 3','Adversary 1','Adversary 2','Adversary 3','Adversary 4','Adversary 5','Adversary 6','Adversary 7','Adversary 8','Adversary 9'};
legend_ann{4,1}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Adversary 1'};
legend_ann{4,3}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Adversary 1','Adversary 2','Adversary 3'};
legend_ann{5,1}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1'};
legend_ann{5,3}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1','Adversary 2','Adversary 3'};
legend_ann{5,4}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1','Adversary 2','Adversary 3','Adversary 4'};
legend_ann{5,9}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1','Adversary 2','Adversary 3','Adversary 4','Adversary 5','Adversary 6','Adversary 7','Adversary 8','Adversary 9'};
legend_ann{5,10}={'Annotator 1','Annotator 2','Annotator 3','Annotator 4','Annotator 5','Adversary 1','Adversary 2','Adversary 3','Adversary 4','Adversary 5','Adversary 6','Adversary 7','Adversary 8','Adversary 9','Adversary 10'};

legend(legend_ann{adversaries(1)-1,max(size(adversaries))})

figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',16,'FontName','Helvetica')

%% to plot the bar graphs with the bar's height corresponding to the annotator's score
figure
number_of_annotators=adversaries(end);
x=repmat(p',1,number_of_annotators);
bar(x,mean_score,'b')
hold on;
z=[zeros(max(size(p)),adversaries(1)-1) mean_score(:,adversaries(1):number_of_annotators)];
bar(x,z,'r')
grid on;
delta=0.0178-0.001*number_of_annotators;
if mod(number_of_annotators,2)==0
pt_distrib=-(number_of_annotators/2)*delta:delta:(number_of_annotators/2)*delta;
i=find(pt_distrib);
adjust=repmat(pt_distrib(i),max(size(p)),1);
else
pt_distrib=-floor(number_of_annotators/2)*delta:delta:floor(number_of_annotators/2)*delta;
adjust=repmat(pt_distrib,max(size(p)),1);
end
errorbar(x+adjust,mean_score,err_bar,'kx')
title('Annotator score variation for various p_a')
xlabel('Adversary p_a - flip probability')
ylabel('Score for every annotator')
legend(legend_ann{adversaries(1)-1,max(size(adversaries))})
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',16,'FontName','Helvetica')