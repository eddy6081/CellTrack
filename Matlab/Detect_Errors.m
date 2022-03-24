%%%%%%%
%Author: Chris Eddy
%Goal: Change cell numbers if necessary following tracking. 
%Function: This creates a text file that keeps track of changes for you
%automatically in the same path that the csv file is found. A new csv file
%with the edited changes is saved in the same path.
%user is prompted to say which cell to change to.
%%%%%%%

disp("The objective of this script is to allow the user, once completed with cell tracking ")
disp("to edit the tracking by detecting overlapping cells, if necessary.")
disp("This will run over all available frames, detecting cells within threshold distance.")

disp("Please select .CSV file that was saved after running cell tracking: ")

[file,path]=uigetfile('*.csv');

T=readtable(fullfile(path,file));

rthresh = 10;

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
%count unique frames
F = unique(T.Frames);

for frame1 = 1:length(F)
    currentframe=F(length(F)-(frame1-1)); %working backwards.
    %find all rows in T that have this frame.
    inds = find(T.Frames==currentframe);
    %for each of these items, we will need their cellnum, xcentroid, and
    %ycentroid
    xs = T.CenterX(inds);%T.CentroidX(inds);
    ys = T.CenterY(inds);%T.CentroidY(inds);
    cnums = T.CellNum(inds);
    %     
    %     %first, lets see if there are any cells with the SAME number on the
    %     %same image.
    %     if length(unique(cnums))<length(cnums)
    %         
    %         %find repeating numbers.
    %         unums=unique(cnums);
    %         Ncount=histc(cnums,unums);
    %         %now, whereever there is Ncount>1 is the same cell number in unums.
    %         av_cells=unums(Ncount>1);
    %         %make sure correct variables exist.
    %         if ~exist('AllFrameStats')
    %             %load the last .mat file
    %             disp('Please load the last .mat file after running PreviousCell_2D.m')
    %             [f,p]=uigetfile('*.mat');
    %             load([p,f])
    %         end
    %         if ~exist('Igray')
    %             disp('Please point to a single image file from the experiment, needed for image sizing')
    %             [f,p]=uigetfile('*.tif');
    %             info=imfinfo([p,f]);
    %             rows=info.Height;
    %             cols=info.Width;
    %         else
    %             [rows,cols]=size(Igray);
    %         end
    %         %now, using the .mat files...
    %         %indsAllA=AllFrames(AllFrameNumbers==av_cells(1));
    %         AS=AllFrameStats(AllFrames==frame1);
    %         statA = AS(
    %         AllStatsA=AllFrameStats(AllFrameNumbers==av_cells(1));
    %         %indsAllB=AllFrames(AllFrameNumbers==av_cells(2));
    %         AllStatsB=AllFrameStats(AllFrameNumbers==av_cells(2));
    %         %shared_frames=indsAllA(ismember(indsAllA,indsAllB));
    %         %if ~isempty(shared_frames)
    %         %    iAllA=ismember(indsAllA,indsAllB);
    %     
    %     %now lets detect if any overlap
    %find if any are within r=10 pixels of each other.
    %create a matrix. 
    %     ind1  ind2  ind3
    %ind1  r=0   
    %ind2        r=0
    %ind3              r=0
    %no, just compute the distance from ind(1) to ind(2:end)
    % then do ind(2) to ind(3:end)
    deleted_list = [];
    for obj = 1:length(inds)-1
        if ~ismember(cnums(obj),deleted_list)
            xdiffs=xs(obj+1:end)-xs(obj);
            ydiffs=ys(obj+1:end)-ys(obj);
            rs = sqrt(xdiffs.^2 + ydiffs.^2);
            %are any of rs less than the threshold? if so, report as possible
            %overlap.
            for i=1:length(rs)
                if rs(i)<rthresh
                    %engage user. 
                    %first, say possible overlap, then bring up tracked images
                    %and plot a large circle around it.
                    %or not, just print the frame and the centroid and have
                    %user go look.
                    disp('Possible overlap detected.')
                    sprintf(['Overlap in frame %d, between cells %d and %d, located at %.2f, %.2f'],currentframe,cnums(obj),cnums(obj+i), xs(obj), ys(obj))
                    av_cells = [cnums(obj),cnums(obj+i)];
                    input_hf1 = str2double(input('Which cell number do you wish to delete? (Enter an integer, or 9999 to combine, 0 to skip)    ', 's'));
                    while isnan(input_hf1) || fix(input_hf1) ~= input_hf1 || ~ismember(input_hf1,[av_cells,0,9999])
                      input_hf1 = str2double(input('Invalid entry. Which cell number do you wish to delete? (Enter an integer, or 9999 to combine, 0 to skip)    ', 's'));
                    end
                    if input_hf1~=0
                        if input_hf1==9999
                            combiner=2;
                        else
                            combiner=1;
                        end
                    else
                        combiner=0;
                    end

                    %if so, change cell number for inds1.
                    if combiner==1
                        T.CellNum(T.CellNum==input_hf1)=av_cells(av_cells~=input_hf1);

                        T=sortrows(T,2,'ascend');

                        fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                        fprintf(fileID,'%d -> %d \n', [input_hf1,av_cells(av_cells~=input_hf1)])
                        fclose(fileID);
                        deleted_list = [deleted_list,av_cells(av_cells~=input_hf1)];
                        disp("Cell number changed.")
                    end

                    if combiner==2
                        %make sure correct variables exist.
                        if ~exist('AllFrameStats')
                            %load the last .mat file
                            disp('Please load the last .mat file after running PreviousCell_2D.m')
                            [f,p]=uigetfile('*.mat');
                            load([p,f])
                        end
                        if ~exist('Igray')
                            disp('Please point to a single image file from the experiment, needed for image sizing')
                            [f,p]=uigetfile('*.tif');
                            info=imfinfo([p,f]);
                            rows=info.Height;
                            cols=info.Width;
                        else
                            [rows,cols]=size(Igray);
                        end
                        %now, using the .mat files...
                        indsAllA=AllFrames(AllFrameNumbers==av_cells(1));
                        AllStatsA=AllFrameStats(AllFrameNumbers==av_cells(1));
                        indsAllB=AllFrames(AllFrameNumbers==av_cells(2));
                        AllStatsB=AllFrameStats(AllFrameNumbers==av_cells(2));
                        shared_frames=indsAllA(ismember(indsAllA,indsAllB));
                        if ~isempty(shared_frames)
                            iAllA=ismember(indsAllA,indsAllB);
                            iAllB=ismember(indsAllB,indsAllA);
                            statsA=AllStatsA(iAllA);
                            statsB=AllStatsB(iAllB);
                            shared_frames_del=[];
                            for frame=1:sum(iAllA)
                                if sum(ismember(statsA(frame).PixelIdxList, statsB(frame).PixelIdxList))>0
                                    bw_temp=zeros(rows,cols);
                                    bw_temp(statsA(frame).PixelIdxList)=1;
                                    bw_temp(statsB(frame).PixelIdxList)=1;
                                    if frame==1 || ~exist('sharestat')
                                        %sharecenters = Max_Fit_Circles(bw_temp);
                                        CC = bwconncomp(bw_temp>0,4);
                                        sharestat=regionprops(CC, {'Area','Centroid',...
                                            'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                                            'PerimeterOld','FilledArea',...
                                            'MajorAxisLength','MinorAxisLength','Eccentricity',...
                                            'Solidity','Extent','ConvexArea','PixelIdxList'});
                                        if length(sharestat)>1
                                            disp('two objects, not just one.')
                                            disp('pause')
                                        end

                                        shareconvperim = Convex_Perimeter(bw_temp>0);
                                        shareFL = Fiber_Length(bw_temp>0);
                                        [sharecenters,shareradii] = Max_Fit_Circles(bw_temp>0); %see function
                                        shareBM = Bleb_Length(bw_temp>0,sharecenters,shareradii);
                                    else
                                        CC = bwconncomp(bw_temp>0,4);
                                        stat=regionprops(CC, {'Area','Centroid',...
                                            'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
                                            'PerimeterOld','FilledArea',...
                                            'MajorAxisLength','MinorAxisLength','Eccentricity',...
                                            'Solidity','Extent','ConvexArea','PixelIdxList'});
                                        if length(stat)>1
                                            disp('Two objects, not just one')
                                            disp('pause')
                                        end

                                        convperim = Convex_Perimeter(bw_temp>0);
                                        FL = Fiber_Length(bw_temp>0);
                                        [centers,radii] = Max_Fit_Circles(bw_temp>0); %see function
                                        BM = Bleb_Length(bw_temp>0,sharecenters,shareradii);
                                        sharestat=[sharestat;stat];
                                        sharecenters = [sharecenters;centers];
                                        shareradii = [shareradii;radii];
                                        shareconvperim = [shareconvperim;convperim];
                                        shareFL = [shareFL;FL];
                                        shareBM = [shareBM;BM];
                                    end
                                else
                                    fprintf('Cells do not have any pixel overlap in frame %d\n', shared_frames(frame))
                                    shared_frames_del(end+1)=frame;
                                end
                            end
                            shared_frames(shared_frames_del)=[];
                            if ~isempty(shared_frames)

                                %save information in CSV File format. 
                                %frame number
                                %correctedNum
                                %area
                                %centroidx
                                %centroidy
                                %major
                                %minor
                                %eccentricity
                                %orientation
                                %convex area
                                %filledarea
                                %equivDiameter
                                %solidity
                                %extent
                                %perimeter
                                %perimeter old
                                %convex perimeter
                                %fiber length
                                %max in radii
                                %bleb length
                                %centersx
                                %centersy

                                A=cat(1,sharestat.Area);
                                cent=cat(1,sharestat.Centroid);
                                posx=cent(:,1);
                                posy=cent(:,2);
                                mjr=cat(1,sharestat.MajorAxisLength);
                                mnr=cat(1,sharestat.MinorAxisLength);
                                ecc=cat(1,sharestat.Eccentricity);
                                ori=cat(1,sharestat.Orientation);
                                cA=cat(1,sharestat.ConvexArea);
                                fA=cat(1,sharestat.FilledArea);
                                eD=cat(1,sharestat.EquivDiameter);
                                s=cat(1,sharestat.Solidity);
                                ex=cat(1,sharestat.Extent);
                                P=cat(1,sharestat.Perimeter);
                                Pold=cat(1,sharestat.PerimeterOld);
                                cP = shareconvperim;
                                fl = shareFL;
                                rad = shareradii;
                                bm = shareBM;
                                cx = sharecenters(:,1);
                                cy = sharecenters(:,2);

                                T2=table(shared_frames,repmat(av_cells(1),length(shared_frames),1), A, posx, posy, mjr, mnr, ecc, ori, ...
                                    cA, fA, eD, s, ex, P, Pold, cP, fl, rad, bm, cx, cy);

                                T2.Properties.VariableNames={'Frames','CellNum','Area','CentroidX',...
                                    'CentroidY','MjrAxisLength','MnrAxisLength','Eccentricity',...
                                    'Orientation','ConvexArea','FilledArea','EquivDiameter',...
                                    'Solidity','Extent','Perimeter','PerimOld', 'ConvexPerim',...
                                    'FibLen', 'InscribeR','BlebLen','CenterX','CenterY'};

                                T2=sortrows(T2,2,'ascend');

                                %from T1 now, delete all elements where the two cells
                                %had overlapping frames.

                                for frame=1:length(shared_frames)
                                    indsA=find(T.Frames==shared_frames(frame) & T.CellNum==av_cells(1));
                                    T(indsA,:)=[];
                                    indsB=find(T.Frames==shared_frames(frame) & T.CellNum==av_cells(2));
                                    T(indsB,:)=[];
                                end

                                T=[T;T2];
                                T=sortrows(T,2,'ascend');

                                fileID=fopen([path,'EditedcsvREADME.txt'],'a');
                                fprintf(fileID,'%d combined with %d \n', [av_cells(2),av_cells(1)])
                                fclose(fileID);
                                deleted_list = [deleted_list,av_cells(2)];
                            else
                                disp('Cells not combined, no overlapping frames')
                            end
                            %now find where these are at in the table. 
                            %so for each of these frames we want to change B->A from
                            %table. Ultimately, we want to combine these pixels and
                            %rerun region props.
                        else
                            disp('These cells do not display any frame overlap')
                        end
                    end
                end
            end
        end
    end
                
        
end



idcs=strfind(file,'.');
writetable(T,[path,file(1:idcs(end)-1),'_edited.csv'],'Delimiter',',','WriteRowNames',true)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm
function [Centers,Rs] = Max_Fit_Circles(BW)
    B=bwboundaries(BW,4,'noholes'); %cell array
    [objs,~]=size(B);
    Centers = zeros(objs,2);
    Rs = zeros(objs,1);
    for bb=1:objs
        grain=false(size(BW));
        for uu=1:length(B{bb})
            grain(B{bb}(uu,1),B{bb}(uu,2))=true;
        end
        grain=uint8(grain)*255;
        [R,cx,cy] = max_inscribed_circle(grain,[]);
        Centers(bb,1)=cx;
        Centers(bb,2)=cy;
        Rs(bb) = R;
    end
end


function ConvP = Convex_Perimeter(BW)
    %the bw image should already be filtered by size.
    %from black and white, create concave objects.
    bw_labeled=bwlabel(BW,4);
    ConvP=zeros(max(max(bw_labeled)),1);
    %we might think that doing bwconvhull(bw,'objects') is indeed faster,
    %but it fails when two objects are nearby and when turned convex they
    %overlap.
    for i=1:max(max(bw_labeled))
        objs = bwconvhull(bw_labeled==i);
        %calculate the perimeter.
        CC = bwconncomp(objs>0,4);
        stats = regionprops(CC, {'Perimeter'});
        ConvP(i)=stats.Perimeter;
    end
end

function FL = Fiber_Length(BW)
    %the bw image should already be filtered by size and in correct binary
    %form.
    bw2=bwlabel(BW,4);
    FL = zeros(max(max(bw2)),1);
    for obj=1:max(max(bw2))
        bw1 = bw2==obj;
        %calculate longest skeleton path, aka fiber length for each object
        bw_skel=bwmorph(bw1,'skel','Inf');
        %how many objects are there?
        ep = find(bwmorph(bw_skel,'endpoints'));
        fl=0;
        for d=1:length(ep)
            D=bwdistgeodesic(bw_skel,ep(d),'quasi');
            if max(max(D))>fl
                fl=max(max(D));
            end
        end
        FL(obj) = fl;
    end
end

function max_bleb_len = Bleb_Length(BW,Centers,Rs)
    bwl = bwlabel(BW,4);
    max_bleb_len = zeros(max(max(bwl)),1);
    for obj=1:max(max(bwl))
        %find skeleton of obj
        bw = bwl==obj;
        bwskel = bwmorph(bw,'skel','Inf');
        Itemp = zeros(size(BW));
        Itemp(Centers(obj,2),Centers(obj,1))=1;
        Id = bwdist(Itemp);
        bw2 = bwskel;
        bw2(Id<ceil(Rs(obj)))=0;
        %find endpoints
        ep = find(bwmorph(bw2,'endpoints'));
        max_d = 0;
        for endpoint = 1:length(ep)
            bw3 = bwdistgeodesic(bw2,ep(endpoint),'quasi');
            m = nanmax(bw3(~isinf(bw3)));
            if m > max_d
                max_d = m; 
            end
        end
        max_bleb_len(obj) = max_d;
    end
    %     if max(max(bwl))==1
    %         C = imfuse(BW,bwskel,'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);
    %         imshow(C)
    %         hold on
    %         viscircles(Centers,Rs);
    %     else
    %         bwskel = bwmorph(BW,'skel','Inf');
    %         C = imfuse(BW,bwskel,'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);
    %         imshow(C)
    %         hold on
    %         viscircles(Centers,Rs);
    %     end
    %could also consider measuring bleb coherence, that is, the
    %directionality of all protrusions - random, in line or one sided.
end
