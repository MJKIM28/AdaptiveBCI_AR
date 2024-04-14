function nasatlx_ButtonDown_ClickRate(hObject, eventdata)
global axesloc
global response_rating
global objmade

    aa = get(gca,'CurrentPoint');
    score = 5*round(aa(1)/5);
    aaa = get(gca);
    axesselected = find(axesloc == aaa.Position(2)); % find which question answered
    
    if ~isempty(response_rating(axesselected)) % if this question already answered 
        delete(objmade{axesselected}); % erase previous answer      
    end

    response_rating(axesselected) = score;  
    objmade{axesselected} = plot(gca,score,0.5,'rv','Markersize',5,'linewidth',5); % update the answer
    

end