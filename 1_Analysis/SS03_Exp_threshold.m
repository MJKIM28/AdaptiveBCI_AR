close all; clear all

%%

LAMBDA1 = 0;
LAMBDA2 = 0;

Thresholds = [0.2:0.2:2];
%%

for threshold1 =Thresholds
    for threshold2 = Thresholds
SS03_main(LAMBDA1,LAMBDA2,threshold1,threshold2)
    end
end
