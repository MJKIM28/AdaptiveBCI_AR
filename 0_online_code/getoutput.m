function output = getoutput(score)
% score: trial x stim
[~,output] = max(sum(score));

end