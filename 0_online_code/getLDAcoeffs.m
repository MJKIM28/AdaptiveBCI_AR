function Coeffs = getLDAcoeffs(mdl)
Diff = -diff(mdl.Mu,1);
w = mdl.invSigma*Diff';
b = -w'*mdl.MuTotal'+log(mdl.N(1)/mdl.N(2));

% nontarget
Coeffs(1).Linear = w;
Coeffs(1).Const = b;

% target
Coeffs(2).Linear = -w;
Coeffs(2).Const = -b;