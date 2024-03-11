function [Decision,class] = decidestopping(pvalue, alpha)


c = getclass(pvalue,alpha);
if c~= -1
    Decision = 1;
    class = c;
else
    Decision = -1;
    class = -1;
end