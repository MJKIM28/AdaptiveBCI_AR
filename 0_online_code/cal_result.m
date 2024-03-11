function cal_result(type)
Ntry = input('# tried:');
Ncorrect = input('# correct:');

if strcmp(type,'AIRD')
    
    acc = Ncorrect/Ntry;
    N = 3; repeat = 10; isi= 0.1*2;
    
elseif strcmp(type,'FAN')
    
    acc = Ncorrect/Ntry;
    N = 6; repeat = 10; isi= 0.1*2;
    
elseif strcmp(type,'Def')
    
    acc = Ncorrect/Ntry;
    N = 4; repeat = 5; isi= 0.075*2;
    
end

if acc == 1
    acc_new = 0.999999;
    ITR = (log2(N) + acc_new*log2(acc_new) + (1-acc_new)*log2((1-acc_new)/(N-1)))*60/(4*repeat*isi);
else
    ITR = (log2(N) + acc*log2(acc) + (1-acc)*log2((1-acc)/(N-1)))*60/(4*repeat*isi);
end
fprintf('Acc %.2f ITR %.2f\n',acc,ITR);