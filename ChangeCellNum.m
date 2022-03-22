%%%%%%%
%Author: Chris Eddy
%Goal: Change cell numbers if necessary following tracking. 
%Function: This creates a text file that keeps track of changes for you
%automatically in the same path that the csv file is found. A new csv file
%with the edited changes is saved in the same path.
%user is prompted to say which cell to change to.
%%%%%%%

disp("The objective of this script is to allow the user, once completed with cell tracking ")
disp("to edit the tracking and combine cells that cell numbers if necessary.")
disp("This will run recursively until the user chooses to exit following the prompt.")

disp("Please select .CSV file that was saved after running cell tracking: ")

[file,path]=uigetfile('*.csv');

T=readtable(fullfile(path,file));

%generate text file to keep track of edited changes. 
if ~exist([path, 'EditedcsvREADME.txt'])
    fileID=fopen([path,'EditedcsvREADME.txt'],'w');
    fprintf(fileID,'This was automatically generated to track edits to original CSV file.\n');
    fclose(fileID);
end

fileID=fopen([path,'EditedcsvREADME.txt'],'a');
fprintf(fileID,'Changes made on %s listed below: \n', string(datetime('now')))
fclose(fileID);

input_hf1=2;
while input_hf1~=0
    %first form unique list.
    cellnums=unique(T.CellNum);
    frames=T.Frames;

    input_hf1 = str2double(input('Which cell number do you wish to change? (Enter an integer, 0 to exit)    ', 's'));
    while isnan(input_hf1) || fix(input_hf1) ~= input_hf1 || ~ismember(input_hf1,[cellnums;0])
      input_hf1 = str2double(input('Invalid entry. Which cell number do you wish to change? (Enter an integer, 0 to exit)    ', 's'));
    end
    
    if input_hf1~=0
        input_hf2 = str2double(input('Which cell number do you wish to change it to? (Enter an integer, 0 to delete)    ', 's'));
        while isnan(input_hf2) || fix(input_hf2) ~= input_hf2 || input_hf2==input_hf1 %|| ~ismember(input_hf2,[cellnums;0;-1])
          input_hf2 = str2double(input('Invalid entry. Which cell number do you wish to change it to? (Enter an integer, 0 to delete, -1 to skip)    ', 's'));
        end
        
        if input_hf2>0
            inds1=find(T.CellNum==input_hf1);
            inds2=find(T.CellNum==input_hf2);

            combiner=1;
            %%%CHECK IF THESE CELLS CAN BE COMBINED.
            if ismember(frames(inds1),frames(inds2))
                a=frames(inds1);
                b=frames(inds2);
                match_ind=ismember(a,b);
                for i=1:length(nonzeros(match_ind))
                    fprintf('Cannot combine these cells, they were simultaneously detected in frame %d. \n', a(i))
                end
                combiner=0;
                
                input_hf3 = str2double(input('Do you wish to combine these cells into a single cell? (Enter an integer) 1: yes, 2: no, 3: skip    ', 's'));
                while isnan(input_hf3) || fix(input_hf3) ~= input_hf3 || input_hf3>3 || input_hf3<1
                  input_hf3 = str2double(input('Invalid entry. Do you wish to combine these cells into a single cell? (Enter an integer) 1: yes, 2: no, 3: skip    ', 's'));
                end
                
                if input_hf3==1
                    combiner=2;
                end
            end
            
            input_hf4 = str2double(input('Do you wish to change only particular frames? (Enter an integer) 1: yes, 2: no (combine all frames), 3: skip (reenter)    ', 's'));
            while isnan(input_hf4) || fix(input_hf4) ~= input_hf4 || input_hf4>3 || input_hf4<1
              input_hf4 = str2double(input('Invalid entry. Do you wish to change only particular frames? (Enter an integer) 1: yes, 2: no (combine all frames), 3: skip (reenter)    ', 's'));
            end
            if input_hf4==1
                combiner=3;
                input_hf5 = input('Please enter a single frame number, or a sequence as [1:3,6,7:9]    ');
                %inds1=find(T.CellNum==input_hf1);
                indstochange1=[];
                for entry = 1:length(inds1)
                    if ismember(T.Frames(inds1(entry)),input_hf5)
                        indstochange1=[indstochange1;inds1(entry)];
                    end
                end
                indstochange2=[];
                for entry = 1:length(inds2)
                    if ismember(T.Frames(inds2(entry)),input_hf5)
                        indstochange2=[indstochange2;inds2(entry)];
                    end
                end
                T.CellNum(indstochange1)=input_hf2;
                T.CellNum(indstochange2)=input_hf1;

                T=sortrows(T,2,'ascend');

                fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                if ~isempty(indstochange1) && ~isempty(indstochange2) %both are not empty.
                    fprintf(fileID,'switched %d <-> %d in frames %s \n', [input_hf1,input_hf2,mat2str(input_hf5)])
                elseif isempty(indstochange2)
                    fprintf(fileID,'switched %d -> %d in frames %s \n', [input_hf1,input_hf2,mat2str(input_hf5)])
                end
                fclose(fileID);
                disp("Cell number changed.")
            
            elseif input_hf4==3
                disp('skipping')
                combiner=5;
            end

            %if so, change cell number for inds1.
            if combiner==1
                T.CellNum(inds1)=input_hf2;

                T=sortrows(T,2,'ascend');
                
                fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                fprintf(fileID,'%d -> %d \n', [input_hf1,input_hf2])
                fclose(fileID);
                disp("Cell number changed.")
            end
            
            if combiner==2
                disp('Need to edit this code, use Detect_Errors.m')
                %                 %make sure correct variables exist.
                %                 if ~exist('AllFrameStats')
                %                     %load the last .mat file
                %                     disp('Please load the last .mat file after running PreviousCell_2D.m')
                %                     [f,p]=uigetfile('*.mat');
                %                     load([p,f])
                %                 end
                %                 if ~exist('Igray')
                %                     disp('Please point to a single image file from the experiment, needed for image sizing')
                %                     [f,p]=uigetfile('*.tif');
                %                     info=imfinfo([p,f]);
                %                     rows=info.Height;
                %                     cols=info.Width;
                %                 else
                %                     [rows,cols]=size(Igray);
                %                 end
                %                 %now, using the .mat files...
                %                 indsAllA=AllFrames(AllFrameNumbers==input_hf1);
                %                 AllStatsA=AllFrameStats(AllFrameNumbers==input_hf1);
                %                 indsAllB=AllFrames(AllFrameNumbers==input_hf2);
                %                 AllStatsB=AllFrameStats(AllFrameNumbers==input_hf2);
                %                 shared_frames=indsAllA(ismember(indsAllA,indsAllB));
                %                 if ~isempty(shared_frames)
                %                     iAllA=ismember(indsAllA,indsAllB);
                %                     iAllB=ismember(indsAllB,indsAllA);
                %                     statsA=AllStatsA(iAllA);
                %                     statsB=AllStatsB(iAllA);
                %                     for frame=1:sum(iAllA)
                %                         bw_temp=zeros(rows,cols);
                %                         bw_temp(statsA(frame).PixelIdxList)=1;
                %                         bw_temp(statsB(frame).PixelIdxList)=1;
                %                         if frame==1
                %                             sharecenters = Max_Fit_Circles(bw_temp);
                %                             sharestat=regionprops(bw_temp>0, {'Area','Centroid',...
                %                                 'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                %                                 'PerimeterOld','FilledArea',...
                %                                 'MajorAxisLength','MinorAxisLength','Eccentricity',...
                %                                 'Solidity','Extent','ConvexArea','PixelIdxList'});
                %                             if length(sharestat)>1
                %                                 disp('two objects, not just one.')
                %                                 disp('pause')
                %                             end
                %                         else
                %                             stat=regionprops(bw_temp>0, {'Area','Centroid',...
                %                                 'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                %                                 'PerimeterOld','FilledArea',...
                %                                 'MajorAxisLength','MinorAxisLength','Eccentricity',...
                %                                 'Solidity','Extent','ConvexArea','PixelIdxList'});
                %                             if length(stat)>1
                %                                 disp('Two objects, not just one')
                %                                 disp('pause')
                %                             end
                %                             sharestat=[sharestat;stat];
                %                             sharecenters = [sharecenters;Max_Fit_Circles(bw_temp)];
                %                         end
                %                     end
                %                     
                %                     A=cat(1,sharestat.Area);
                %                     cent=cat(1,sharestat.Centroid);
                %                     posx=cent(:,1);
                %                     posy=cent(:,2);
                %                     mjr=cat(1,sharestat.MajorAxisLength);
                %                     mnr=cat(1,sharestat.MinorAxisLength);
                %                     ecc=cat(1,sharestat.Eccentricity);
                %                     ori=cat(1,sharestat.Orientation);
                %                     cA=cat(1,sharestat.ConvexArea);
                %                     fA=cat(1,sharestat.FilledArea);
                %                     eD=cat(1,sharestat.EquivDiameter);
                %                     s=cat(1,sharestat.Solidity);
                %                     ex=cat(1,sharestat.Extent);
                %                     P=cat(1,sharestat.Perimeter);
                %                     Pold=cat(1,sharestat.PerimeterOld);
                %                     cx = sharecenters(:,1);
                %                     cy = sharecenters(:,2);
                % 
                % 
                %                     T2=table(shared_frames,repmat(input_hf1,length(shared_frames),1), A, posx, posy, mjr, mnr, ecc, ori, ...
                %                         cA, fA, eD, s, ex, P, Pold, cx, cy);
                % 
                %                     T2.Properties.VariableNames={'Frames','CellNum','Area','CentroidX',...
                %                         'CentroidY','MjrAxisLength','MnrAxisLength','Eccentricity',...
                %                         'Orientation','ConvexArea','FilledArea','EquivDiameter',...
                %                         'Solidity','Extent','Pertimeter','PerimOld',...
                %                         'CenterX','CenterY'};
                % 
                %                     T2=sortrows(T2,2,'ascend');
                % 
                %                     %from T1 now, delete all elements where the two cells
                %                     %had overlapping frames.
                %                     
                %                     for frame=1:length(shared_frames)
                %                         indsA=find(T.Frames==shared_frames(frame) & T.CellNum==input_hf1);
                %                         T(indsA,:)=[];
                %                         indsB=find(T.Frames==shared_frames(frame) & T.CellNum==input_hf2);
                %                         T(indsB,:)=[];
                %                     end
                %                     
                %                     T=[T;T2];
                %                     T=sortrows(T,2,'ascend');
                %                     
                %                     fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                %                     fprintf(fileID,'%d combined with %d \n', [input_hf2,input_hf1])
                %                     fclose(fileID);
                %                     %now find where these are at in the table. 
                %                     %so for each of these frames we want to change B->A from
                %                     %table. Ultimately, we want to combine these pixels and
                %                     %rerun region props.
                %                 else
                %                     disp('These cells do not display any frame overlap')
                %                 end
            end
                
        elseif input_hf2==0
            input_hf3 = str2double(input('Do you wish to delete all frames with this cell? (Enter an integer) 1: yes, 2: no    ', 's'));
            while isnan(input_hf3) || fix(input_hf3) ~= input_hf3 || input_hf3>2 || input_hf3<1
              input_hf3 = str2double(input('Invalid entry. Do you wish to delete all frames with this cell? (Enter an integer) 1: yes, 2: no      ', 's'));
            end
            if input_hf3==1
                inds1=find(T.CellNum==input_hf1);
                T(inds1,:)=[];
                fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                fprintf(fileID,'Cell %d deleted \n', input_hf1)
                fclose(fileID);
                disp("Cell number deleted.")
            else
                input_hf4 = input('Please enter a single frame number, or a sequence as [1:3,6,7:9]    ');
                inds1=find(T.CellNum==input_hf1);
                indstodelete=[];
                for entry = 1:length(inds1)
                    if ismember(T.Frames(inds1(entry)),input_hf4)
                        indstodelete=[indstodelete;inds1(entry)];
                    end
                end
                T(indstodelete,:)=[];
                fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                fprintf(fileID,'Cell %d deleted in frame(s) %s \n', [input_hf1,mat2str(input_hf4)])
                fclose(fileID);
                disp("Cell number deleted in selected frames.")
                
            end
        end
    else
        disp("Exiting and saving edited .CSV file in same path...")
    end
end

idcs=strfind(file,'.');
writetable(T,[path,file(1:idcs(end)-1),'_edited.csv'],'Delimiter',',','WriteRowNames',true)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm

function Centers = Max_Fit_Circles(BW)
    B=bwboundaries(BW,'noholes'); %cell array
    [objs,~]=size(B);
    Centers = zeros(objs,2);
    for bb=1:objs
        grain=false(size(BW));
        for uu=1:length(B{bb})
            grain(B{bb}(uu,1),B{bb}(uu,2))=true;
        end
        grain=uint8(grain)*255;
        [~,cx,cy] = max_inscribed_circle(grain,[]);
        Centers(bb,1)=cx;
        Centers(bb,2)=cy;
    end
end