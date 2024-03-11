function sc = predict_mj(mdl,x)
sc = (sum(mysigmoid(mdl.SV,x).*mdl.Y_SV.*mdl.Alpha)+mdl.Bias)';
sc = [-sc sc];
end

