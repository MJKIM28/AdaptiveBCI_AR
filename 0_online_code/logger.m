function logger(str, sess,path)
if nargin < 2
    sess            = 0;
end
dateS = string(datetime('today'));
fileID              = fopen(sprintf([path,'/DB_log/EXP_LOG_%s.txt'],dateS),'a+');
clk                 = string(datetime('now'));
fprintf(fileID,'\r\n[%s] %s: %s',sess, clk, str);
fclose(fileID);