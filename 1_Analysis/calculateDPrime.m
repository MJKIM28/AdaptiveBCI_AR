function d_prime = calculateDPrime(hitRate, falseAlarmRate)
    % Convert rates to Z scores
    zHit = norminv(hitRate);
    zFalseAlarm = norminv(falseAlarmRate);
    
    % Calculate d'
    d_prime = zHit - zFalseAlarm;
end