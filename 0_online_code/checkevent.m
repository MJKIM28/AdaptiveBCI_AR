function result = checkevent(trigger)

% Check trigger event
% whether all trigger types (e.g. 1,2,3,4) appeared

Count = accumarray(trigger(trigger~=0)',1);

if length(unique(Count)) == 1
    result = true;
else
    result = false;
end