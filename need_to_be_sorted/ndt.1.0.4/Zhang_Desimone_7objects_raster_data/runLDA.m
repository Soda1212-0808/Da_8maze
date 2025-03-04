%pull_psth('20170329ZJF#J024OFC4720');

clear;
clc;

pull_psth('20170329ZJF#J024OFC4720');
psth1 = binned_data{1,1};
s = size(psth1);
length_of_psth = s(2); 

label = binned_labels.unpoke;
label = label{1,1};
label = label';

accuracy = zeros(1,s(2));
for i = 1:s(2)
    
    for j = 1:26
        unit = binned_data{1,j};
        spikes(:,j) = unit(:,i);     
    end
    
    [coeff,score,latent] = pca(spikes);
    %score = score*coeff';
    spikes = score(:,1:3);
    
    
    s = size(spikes);
    correct = zeros(1,s(1));
    for k = 1:s(1)
        index = true(s(1),1);
        index(k) = false;
        training_data = spikes(index,:);
        training_label = label(index);
        test_label = label(~index);
        test = spikes(~index,:);
        assignin('base','training_data',training_data);
        assignin('base','training_label',training_label);
        assignin('base','test',test);
        MdlLinear = fitcdiscr(training_data,training_label);
        p1 = predict(MdlLinear,test);
        if strcmp(p1,test_label)
            correct(k) = 1;
        end
    end
    accuracy(i) = sum(correct)*100/length(correct);
end

plot(-0.2:0.1:0.5,accuracy);

