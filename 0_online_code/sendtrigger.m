function sendtrigger(portobj,val)
fwrite(portobj,val,'sync');
pause(0.01);