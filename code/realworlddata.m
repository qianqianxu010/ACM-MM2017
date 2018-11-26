clc;
clear;


input_filename = 'data1.txt';
%data1.txt is ref1 in PC-VQA dataset
%data3_online.txt is ref3 in PC-IQA dataset

data_ref = importdata(input_filename);

options.intercept = 1;
options.alpha = 0.75;  
options.beta1 = 0.8;   
options.beta2 = 1.03;  

[score, output] = AODHodgerank(data_ref, options);
out=output.outlier_detect(:,end);
outlier_detected=[];
for i=1:length(out)
    if out(i)~=0
outlier_detected=[outlier_detected
    data_ref(i,:)];
    end
end
 