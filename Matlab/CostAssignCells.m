function [newcellnums]=CostAssignCells(thisstats, prevstats, prevcellnums, maxcellnumber, cell_size)
    %new method finds distance from all previous cells. 
    %FUNCTION IS CALLED IN CellImageProcessing.m line 278.
    %%%INPUTS%%%
    %s is a list of the number of objects/cells in frame
    %T2 is matrix of table (we will write) from the previous frame
    %stats = regionprops ran on current frame image
    %blackout = vector containing cell numbers that have divided
    %%%OUTPUTS%%%
    %Newlist is a list of corrected cell numbers for given frame
    %blackoutlist is a list of cells to ignore in future frames following
    %division events
    %%%AUTHOR: Chris Eddy%%%%
    %%%DATE: 08/02/22%%%
    disp("CALCULATING COSTS!")
    
    thisdata = struct2cell(thisstats);
    % thisdata(12,:)=[]; %delete pixel index list, unless we do IoU?
    prevdata = struct2cell(prevstats);
    
    costarray = zeros(size(thisdata,2),size(prevdata,2));
    
    for i=1:size(thisdata,2)
        for j = 1:size(prevdata,2)
            %want to compute cost matrix. If the previous frame has P cells, and
            %this frame has N cells, then we should have an N x P matrix. 
            %we will compute the cost as abs(1 - old/(new+eps))?
            for k = [1,3,4,5,7,8,9,10,11,13,16,17,18,19]
                costarray(i,j) = costarray(i,j) + abs(1 - (thisdata{k,i}+eps)/(prevdata{k,j}+eps));
            end
            %centroid only
            costarray(i,j) = costarray(i,j) + sqrt(sum((thisdata{2,i} - prevdata{2,j}).^2))/cell_size; %56 is approximate cell size in pixels for 20x.
            %calculate IoU
            [~,mi] = intersect(thisdata{12,i},prevdata{12,j});
            [~,ui] = union(thisdata{12,i}, prevdata{12,j});
            costarray(i,j) = costarray(i,j) + 1 - size(mi,1)/size(ui,1);
        end
    end
    
    %now run linear sum assignment
    matches = matchpairs(costarray, 1e4); %play around with the second value, it is essentially a hard match requirement at 1e4
    %matches will be N x 2
    disp("DONE MATCHING")
    
    %%%Assign new cell numbers. 
    %prevcellnums should be P x 1
    newcellnums = zeros(size(thisdata,2),1);
    for i=1:size(newcellnums,1)
        if any(matches(:,1)==i)
            %if it is in the list, assign it to the previous cell label
            newcellnums(i) = prevcellnums(matches(matches(:,1)==i,2));
        else
            %make it a new cell to follow.
            newcellnums(i) = maxcellnumber+1;
            maxcellnumber = maxcellnumber+1;
        end
    end
    disp("DONE MAKING CELL NUMBERS!")
    
    %{
    Just a few notes:
    this is not set up to deal with division events well. Basically it is going
    to try and assign one of the two daughter cells to the parantal cell, and
    the other will be marked as a new cell. 
    %}
end