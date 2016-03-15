clear all; close all; clc

set(0, 'DefaultAxesFontSize', 16);

%...setup working directories
addpath ../MultiLabelerMethods/
addpath ../ErrorStatistics/
addpath ../DATA

% Set the desired 'k' (the number of adversaries) here.
% Number of extra annotators, adversaries, introduced
k=1;

% loading Breast Dataset
load breastdata
Data=norm_Data.X;
SampleSize=size(Data,1); FeatureDim=size(Data,2);
Data=Data-ones(SampleSize,1)*mean(Data);
Data=Data./(ones(SampleSize,1)*std(Data));
Anno=[norm_Data.Y' repmat(norm_Data.Y_golden',1,k)]; LabelerNo=size(Anno, 2);
Labels=norm_Data.Y_golden';
clear norm_Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % loading Atrial Fibrillation Dataset
% load AtrialFib_Clean
% Data=Sample;
% SampleSize=size(Data,1); FeatureDim=size(Data,2);
% Data=Data-ones(SampleSize,1)*mean(Data);
% Data=Data./(ones(SampleSize,1)*std(Data));
% Anno=[Doctor repmat(Truth,1,k)]; LabelerNo=size(Anno, 2);
% Labels=Truth;
% clear Sample Doctor Truth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % loading Glass Dataset
% load DATA/multilabelerGlass
% Data=X;
% SampleSize=size(Data,1); FeatureDim=size(Data,2);
% Data=Data-ones(SampleSize,1)*mean(Data);
% Data=Data./(ones(SampleSize,1)*std(Data));
% Anno=[MLabel repmat(originalZ,1,k)]; LabelerNo=size(Anno, 2);
% Labels=originalZ;
% clear X MLabel originalZ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % loading Iono Dataset
% load ionodata_clean
% Data=Sample;
% Data=Data(:,[1 3:end]);
% SampleSize=size(Data,1); FeatureDim=size(Data,2);
% Data=Data-ones(SampleSize,1)*mean(Data);
% Data=Data./(ones(SampleSize,1)*std(Data));
% Anno=[Doctor repmat(Truth,1,k)]; LabelerNo=size(Anno, 2);
% Labels=Truth;
% clear Sample Doctor Truth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % loading Housing Dataset
% load housingdata_multi
% Data=A;
% SampleSize=size(Data,1); FeatureDim=size(Data,2);
% Data=Data-ones(SampleSize,1)*mean(Data);
% Data=Data./(ones(SampleSize,1)*std(Data));
% Anno=[[d{1} d{2} d{3} d{4} d{5}] repmat(gold_d,1,k)]; LabelerNo=size(Anno, 2);
% Labels=gold_d;
% clear A d gold_d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:SampleSize
    for j=1:1:LabelerNo
        if Anno(i,j)==-1
            Anno(i,j)=0;
        end
    end
    if Labels(i,1)==-1
        Labels(i,1)=0;
    end
end

PrPos=size(find(Labels==1),1)/SampleSize;
PrNeg=1-PrPos;
Blocks=5;
ProbLabelGau=zeros(SampleSize,1);
ProbLabelBino=zeros(SampleSize,1);

%% ....test part....%
Portion=floor(SampleSize/Blocks);
%flip probabilities - amount of noise
p=[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
for flip_noise_level=1:1:size(p,2)
    for Ptr=1:1:Blocks
        if Ptr~=Blocks
            TestData=Data((Ptr-1)*Portion+1:1:Ptr*Portion,:);
            TrainData=Data(setdiff((1:1:SampleSize),((Ptr-1)*Portion+1:1:Ptr*Portion)),:);
            TrainAnno=Anno(setdiff((1:1:SampleSize),((Ptr-1)*Portion+1:1:Ptr*Portion)),:);
            
        else
            TestData=Data((Ptr-1)*Portion+1:1:SampleSize,:);
            TrainData=Data(setdiff((1:1:SampleSize),((Ptr-1)*Portion+1:1:SampleSize)),:);
            TrainAnno=Anno(setdiff((1:1:SampleSize),((Ptr-1)*Portion+1:1:SampleSize)),:);
        end
        
        
        % code for flipping labels
        %jth annotator is an adversary
        for adv=(LabelerNo-k+1):LabelerNo
            for i = 1:size(TrainAnno,1)
                %flipping probability
                if binornd(1,p(flip_noise_level)) == 1
                    TrainAnno(i,adv) = xor(TrainAnno(i,adv),1);
                end
            end
        end
        
        
        %...choose gaussian
        [alpha_g, W_g, PtildZ_g, Pz_x_g, maxLogScr_g, Py_xz_g, Py_xzbar_g]=MLLabeler(TrainData, TrainAnno, 1); %...initialize by default
        [alpha_g_rnd, W_g_rnd, PtildZ_g_rnd, Pz_x_g_rnd, maxLogScr_g_rnd]=MLLabeler(TrainData, TrainAnno, 1, 1); %...initialize by random
        if maxLogScr_g_rnd>maxLogScr_g
            alpha_g=alpha_g_rnd;
        end
        
        %...choose binomial
        [alpha_b, W_b, PtildZ_b, Pz_x_b, maxLogScr_b, Py_xz_b, Py_xzbar_b]=MLLabeler(TrainData, TrainAnno, 2); %...initialize by default
        [alpha_b_rnd, W_b_rnd, PtildZ_b_rnd, Pz_x_b_rnd, maxLogScr_b_rnd, Py_xz_b_rnd, Py_xzbar_b_rnd]=MLLabeler(TrainData, TrainAnno, 2, 1); %...initialize by random
        if maxLogScr_b_rnd>maxLogScr_b
            alpha_b=alpha_b_rnd;
        end
        
        %% p_ann computed for each annotator over a fold of the training data. 5 folds in total
        %binomial case, random initialization
        [p_ann_b_rnd,score1_rnd]=adversarial(Pz_x_b_rnd,Py_xzbar_b_rnd,Py_xz_b_rnd,LabelerNo);
        score1_across_folds(Ptr,:)=score1_rnd;
        
        if Ptr~=Blocks
            ProbLabelGau((Ptr-1)*Portion+1:1:Ptr*Portion,1)=MLLabeler_testing(TestData,alpha_g);
            ProbLabelBino((Ptr-1)*Portion+1:1:Ptr*Portion,1)=MLLabeler_testing(TestData,alpha_b);
        else
            ProbLabelGau((Ptr-1)*Portion+1:1:SampleSize,1)=MLLabeler_testing(TestData,alpha_g);
            ProbLabelBino((Ptr-1)*Portion+1:1:SampleSize,1)=MLLabeler_testing(TestData,alpha_b);
        end
    end
    
    AccuracyGau=0;
    AccuracyBino=0;
    for i=1:1:SampleSize
        if Labels(i,1)==1
            if (ProbLabelGau(i,1)>=PrNeg)
                AccuracyGau=AccuracyGau+1;
            end
            if (ProbLabelBino(i,1)>=PrNeg)
                AccuracyBino=AccuracyBino+1;
            end
        else
            if (ProbLabelGau(i,1)<PrNeg)
                AccuracyGau=AccuracyGau+1;
            end
            if (ProbLabelBino(i,1)<PrNeg)
                AccuracyBino=AccuracyBino+1;
            end
        end
    end
    AccuracyGau=AccuracyGau/SampleSize;
    AccuracyBino=AccuracyBino/SampleSize;
    
    [~, detAxeGau]=rocPlot(Labels, ProbLabelGau, 0, 0, 1);
    [faAxe, detAxeBino]=rocPlot(Labels, ProbLabelBino, 0, 0, 1);
    figure(1), plot(faAxe,detAxeGau, '-r'), grid on, hold on;
    figure(1), plot(faAxe,detAxeBino, '--b');
    
    %mean score vs flip prob. with std. error bars for each annotator
    mean_score1(flip_noise_level,:)=mean(score1_across_folds,1);
    error_bar_for_score1(flip_noise_level,:)=std(score1_across_folds,1,1);
    
    auc(:,flip_noise_level)=trapz(faAxe,detAxeBino);
    
end
%plot mean score vs flip prob. with std. error bars for each annotator
%see function definition for understanding the input arguments to the function
plot_scores_vs_p_with_error_bar(mean_score1,error_bar_for_score1,(LabelerNo-k+1):LabelerNo,p)

figure
plot(p,auc,'k+','LineWidth',10)
hold on
plot(p,auc,'r','LineWidth',2)  %plot the area under the roc curve
title('Performance of the Model when Adversarial Noise Added to Ground Truth','FontSize',25,'FontWeight','bold')
xlabel('Amount of Adversarial Noise','FontSize',25,'FontWeight','bold');
ylabel('Area under the ROC curve','FontSize',25,'FontWeight','bold');
axis([0 1 0 1])