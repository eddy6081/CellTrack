%%%
%Output:
%Input:
%Functions: Max_Fit_Circles and PreivousCell3D defined later in code.
%%% PURPOSE: Use in combination with deep learning cellpose output (saved
%%% .mat files) in order to quickly quantify cell shape, feed into svm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_scale=1; %lets just fix this in the SVM script.
%Check that external functions exist within the directory.
DefaultDir=pwd;
if ~exist([DefaultDir, '/max_inscribed_circle.m'])
    error('The Matlab script "max_inscribed_circle.m" is not in the same directory')
end
if ~exist([DefaultDir, '/inpoly.m'])
    error('The Matlab script "inpoly.m" is not in the same directory')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select where you want to save the printed images to.
disp('Select directory you wish to save output .mat files to (csv directory up)')
save_path=uigetdir(pwd,'Select directory you wish to save output .mat files to');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we need user to select where the saved .mat files are.
disp('Select directory where saved .mat containing labels are from cell_v2.py')
mat_path=uigetdir(pwd,'Select directory where saved .mat containing labels are');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Code to process images and correctly identify cell numbers. 

%for each path, we want to analyze the ith frame, measuring each object and
%assigning a z along with it. These objects have cell numbers (indices)
%assigned to them depending only on where they appear on the image (l->r) by default.
%We want to correct these numbers so that we can accurately track the
%trajectory of any individual cell. 

d_ims=dir([mat_path,filesep,'*.mat']);
d_ims=d_ims(~ismember({d_ims.name},{'.','..'})); %delete hidden files
Matnames=d_ims;

cellmax=0;
maxframe=length(Matnames); %should double check all these filepaths have the same number of files in them

blackoutlist=[]; %%%%%%%%%%%VERY IMPORTANT
allcellnums=[];

%ask user if they want to begin at a certain frame. 
input_hf = str2double(input('Do you wish to begin at a certain frame? (Enter an integer) 1: yes, 2: no    ', 's'));
while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 || input_hf>2
  input_hf = str2double(input('Do you wish to begin at a certain frame? (Enter an integer) 1: yes, 2: no    ', 's'));
end
if input_hf==1
    %uigetfile
    disp('Please guide to the previously saved .mat file')
    [file,path]=uigetfile('*.mat');
    load(fullfile(path,file))
    %all filenames should follow workspace_checkpoint_#.mat so use idcs to
    %look for _ and .mat
    idcPeriod=strfind(file,'.');
    idcUnder=strfind(file,'_');
    AllFrames = cat(1,AllFrameStats.Frame);
    startframe=max(AllFrames)+1;
    %str2double(file(idcUnder(end)+1:idcPeriod(1)-1));
    prev_stats=AllFrameStats(AllFrames==startframe-1);
    prev_cellnums=AllFrameNumbers(AllFrames==startframe-1);
    allcellnums=AllFrameNumbers(AllFrames==startframe-1);
else
    startframe=1;
end

%ask user if they want to begin at a certain frame. 
check_bright_cells = str2double(input('Do you wish to screen cell manually - new cells will be asked to be confirmed by user? (Enter an integer) 1: yes, 2: no    ', 's'));
while isnan(check_bright_cells) || fix(check_bright_cells) ~= check_bright_cells || check_bright_cells<1 || check_bright_cells>2
  check_bright_cells = str2double(input('Invalid entry. Do you wish to screen cell manually - new cells will be asked to be confirmed by user?    ', 's'));
end

for im_num=startframe:maxframe
    %load BW, make logical, run stats, positions, circles, frame number,
    %get cell numbers (add cellmax), assign z
    fprintf('Reading %s\n',Matnames(im_num).name)
    load([Matnames(im_num).folder,filesep,Matnames(im_num).name])
    BW=labels;%(100:2700,300:2800); %labels is actually loaded as an RGB
    
    %by cutting, it would be best to relabel.
    %%%We do this in case the labels are not consecutive. Notice we do not
    %%%relabel the background (0) just the foreground labels. 
    good_labs = unique(BW);
    relabels = [0;good_labs(2:end)-(cumsum(diff(good_labs)-1))];
    for old_i = 2:length(good_labs)
        BW(BW==good_labs(old_i))=relabels(old_i);
    end
    
    %finally, it would be good to remove any pixels that are not connected.
    %the reason is, regionprops perimeter function returns unusual results
    %for discontiguous objects.
    %not computationally efficient, but do bwareaopen for each object.
    good_labs=unique(BW);
    rmd = 0;
    ttl = length(good_labs)-1;
    for old_i=2:length(good_labs)
        F = BW==good_labs(old_i);
        obj_stats = regionprops(F, {'Area', 'PixelIdxList'});
        areas=[obj_stats.Area];
        if max(areas)<50
            for id=1:length(areas)
                BW(obj_stats(id).PixelIdxList)=0;
                rmd=rmd+1;
            end
        else
            if size(obj_stats,1)>1
                %find which has the largest area.
                [~,idx]=max(areas);
                for id=1:length(areas)
                    if id~=idx
                        BW(obj_stats(id).PixelIdxList)=0;
                        rmd=rmd+1;
                    end
                end
            end
        end
    end
    fprintf("removed %d of %d objects.\n",rmd,ttl)
    
    %by cutting, it would be best to relabel.
    %%That is, relabel again in case we removed any detections too small
    good_labs = unique(BW);
    relabels = [0;good_labs(2:end)-(cumsum(diff(good_labs)-1))];
    for old_i = 2:length(good_labs)
        BW(BW==good_labs(old_i))=relabels(old_i);
    end
    
    [h,w,ch]=size(BW);
    %maybe remove objects that are too big?
    %BW=BW-bwareaopen(BW,h*w/80);
    
    %relabel once more. Every channel needs to begin at 1.
    %The current output has only one cell per channel, but still, 
    %at this point we basically make the images binary. However, it still 
    %allows multiple detections to be in a single channel, so more
    %adaptable in case code changes in the future. 
    for c=1:ch
        good_labs = unique(BW(:,:,c));
        relabels = [0;good_labs(2:end)-(cumsum(diff(good_labs)-1))];
        for old_i = 2:length(good_labs)
            BW(BW==good_labs(old_i))=relabels(old_i);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%NOW DO MEASURING%%%%%%%%%%%%%%%
    for c=1:ch
        stats = regionprops(BW(:,:,c), {'Area','Centroid',...
            'Orientation','FilledArea','EquivDiameter', 'Perimeter'...
            'PerimeterOld','FilledArea',...
            'MajorAxisLength','MinorAxisLength','Eccentricity',...
            'Solidity','Extent','ConvexArea','PixelIdxList'});
        if c>1
            stats=[prevstat;stats];
        end
        prevstat = stats;
    end
    pos = cat(1,stats.Centroid);
    %     z = repmat(im_num,size(pos,1),1);
    %     if im_num>1
    %         cellnums=[1:length(z)]'+max(prev_cellnums);
    %     else
    %         cellnums=[1:length(z)]';
    %     end

    %calculate convex hull for convex perimeter
    for c=1:ch
        convperim = Convex_Perimeter(BW(:,:,c));

        [centers,radii] = Max_Fit_Circles(BW(:,:,c)); %see function

        FL = Fiber_Length(BW(:,:,c));

        BM = Bleb_Length(BW(:,:,c),centers,radii);
        
        if c>1
            %Do concatenations
            convperim = [prevconv;convperim];
            centers = [prevcenters;centers];
            radii = [prevradii;radii];
            FL = [prevFL;FL];
            BM = [prevBM;BM];
        end
        prevconv = convperim;
        prevcenters = centers;
        prevradii = radii;
        prevFL = FL;
        prevBM = BM;
    end
            

    OneFrameStats=stats;
    OneFramePos=pos;
    %OneFrameCellnums=cellnums;
    OneFrameCenters=centers;
    OneFrameRadii = radii;
    OneFrameFL = FL;
    OneFrameBM = BM;
    OneFrameConvPerims=convperim;
    
    %throw it all into the stats structure. 
    for i=1:size(stats,1)
        OneFrameStats(i).MICCenter = OneFrameCenters(i,:);
        OneFrameStats(i).MICRadii = OneFrameRadii(i);
        OneFrameStats(i).FiberLength = OneFrameFL(i);
        OneFrameStats(i).BlebLength = OneFrameBM(i);
        OneFrameStats(i).ConvexPerim = OneFrameConvPerims(i);
    end
    
    %%%Now do tracking.
    if im_num>1
        maxcellnum=max(AllFrameNumbers);
        if ~isempty(prev_cellnums)
            [OneFrameCellnums] = CostAssignCells(OneFrameStats, prev_stats, prev_cellnums, maxcellnum, 56);
            l=unique(OneFrameCellnums);
            for cnum=1:length(l)
                if sum(OneFrameCellnums==l(cnum))>1
                    disp('two or more cells display the same corrected cell number')
                end
            end
        else
            OneFrameCellnums=[1:size(OneFrameStats,1)]'+maxcellnum;
        end
    else
        OneFrameCellnums = [1:size(OneFrameStats,1)]';
    end

    %%LOAD IMAGE
    Igray=image;

    %if there are new cells, we want to make sure the user identifies those
    %they do not want to follow.
    %how do you test if there are new cell numbers
    allcellnumspost=unique([OneFrameCellnums;allcellnums]);
    
    if check_bright_cells==1
        if length(allcellnumspost)~=length(allcellnums)
            %we have variables Igray, OneFrameCellnums, OneFrameStats, etc and
            %prev_cellnums, allcellnums,

            BW=zeros(size(Igray));
            T=zeros(size(Igray));
            blackoutpixels=zeros(size(Igray));
            bT=zeros(size(Igray));

            for obj = 1:length(OneFrameCellnums)
                if im_num>1
                    if sum(prev_cellnums==OneFrameCellnums(obj))==0
                        %that is if it is a new cell.
                        BW(OneFrameStats(obj).PixelIdxList)=1;
                        T(OneFrameStats(obj).PixelIdxList)=0.3;
                    end
                else
                    BW(OneFrameStats(obj).PixelIdxList)=1;
                    T(OneFrameStats(obj).PixelIdxList)=0.3;
                end
            end

            x=0;
            y=0;

            while(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                disp('Please select green colored cells you wish to not be followed. If you may click the red cell to change it back to green. Green cells will be tracked.')
                figure(3000)
                imshow(Igray)
                hold on
                gdisp=imshow(cat(3,zeros(size(BW)),BW,zeros(size(BW))));
                set(gdisp,'AlphaData',T)
                hold on
                rdisp=imshow(cat(3,blackoutpixels,zeros(size(blackoutpixels)),zeros(size(blackoutpixels))));
                set(rdisp,'AlphaData',bT)
                title(sprintf('%s',Matnames(im_num).name(1:end-4)))

                [x,y]=ginput(1);
                x=round(x);
                y=round(y);
                if(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                    lindex=sub2ind(size(BW),y,x); %form linear index

                    if BW(lindex)==1 || blackoutpixels(lindex)==1 %that is, if the selected point actually matters
                        %determine which cell number it is.
                        if im_num>1
                            possible_cellnums=setdiff(OneFrameCellnums,prev_cellnums); %this returns the values 
                        else
                            possible_cellnums=OneFrameCellnums;
                        end
                        %in OneFrameCellnums that are not in prev_cellnums
                        %determine which one of these cells were selected.
                        %inds is the cell number.
                        for i=1:length(possible_cellnums)
                            in=find(OneFrameCellnums==possible_cellnums(i)); %should only be length 1
                            if ~isempty(find(OneFrameStats(in).PixelIdxList==lindex))
                                %if the selected point is inside that cell
                                %add to blackoutlist, unless already on it. 
                                if sum(blackoutlist==possible_cellnums(i))<1
                                    %it is not in the list yet.
                                    blackoutlist=[blackoutlist,possible_cellnums(i)];
                                    BW(OneFrameStats(in).PixelIdxList)=0;
                                    %Turn display to red.
                                    T(OneFrameStats(in).PixelIdxList)=0;
                                    blackoutpixels(OneFrameStats(in).PixelIdxList)=1;
                                    bT(OneFrameStats(in).PixelIdxList)=0.3;
                                else
                                    %it is already in the list, delete it.
                                    blackoutlist(blackoutlist==possible_cellnums(i))=[];
                                    %Turn display back to green
                                    BW(OneFrameStats(in).PixelIdxList)=1;
                                    T(OneFrameStats(in).PixelIdxList)=0.3;
                                    blackoutpixels(OneFrameStats(in).PixelIdxList)=0;
                                    bT(OneFrameStats(in).PixelIdxList)=0;
                                end

                                break
                            end
                        end
                    else
                        if(x>-1 && x<size(BW,2)+1 && y>-1 && y<size(BW,1)+1)
                            disp('Invalid point.')
                        end
                    end
                end
            end
        end
    end
    
    allcellnums=OneFrameCellnums;
    
    %store for later
    prev_stats=OneFrameStats;
    prev_cellnums=OneFrameCellnums;
    
    figure
    %Igray=imread([gray_path,filesep,gray_ims(im_num).name]);
    imshow(Igray)
    if ~isempty(OneFrameCellnums)
        for k=1:size(OneFramePos,1)
            if sum(blackoutlist==OneFrameCellnums(k))==0
                c=OneFrameCenters(k,:);
                text(c(1),c(2),sprintf('%d',OneFrameCellnums(k)),...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','middle',...
                    'FontWeight','bold',...
                    'Color','r');
            end
        end
    end
    axis off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=getframe;
    [framerows,framecols,~]=size(h.cdata);

    %%%FIXING OUTPUT SIZE ISSUE.
    if framerows==721
        im=h.cdata(2:end,:,:);
    else
        im=h.cdata;
    end
    %%%

    close all

    imwrite(im,[save_path,filesep,strcat(Matnames(im_num).name(1:end-4),'.tif')])    

    
    frames = repmat(im_num,size(OneFrameStats,1),1);
    for i=1:size(OneFrameStats,1)
        OneFrameStats(i).Frame = frames(i);
    end
    %store all information for the current frame into larger complex
    if im_num==1
        AllFrameStats=OneFrameStats;
        AllFrameNumbers = OneFrameCellnums;
    else
        AllFrameStats=[AllFrameStats; OneFrameStats];
        AllFrameNumbers = [AllFrameNumbers; OneFrameCellnums];
    end

    if mod(im_num,10)==0 || im_num==1
        save(sprintf([save_path,filesep,'workspace_checkpoint_%d.mat'],im_num),'AllFrameStats','AllFrameNumbers','blackoutlist')
    end

end

%Delete the blackout cells we didn't want to track any more.
for obj=1:length(blackoutlist)
    AllFrameStats(AllFrameNumbers==blackoutlist(obj))=[];
    AllFrameNumbers(AllFrameNumbers==blackoutlist(obj))=[];
end

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

A=cat(1,AllFrameStats.Area) * size_scale * size_scale;
Pos = cat(1,AllFrameStats.Centroid);
posx=Pos(:,1);
posy=Pos(:,2);
mjr=cat(1,AllFrameStats.MajorAxisLength) * size_scale;
mnr=cat(1,AllFrameStats.MinorAxisLength) * size_scale;
ecc=cat(1,AllFrameStats.Eccentricity);
ori=cat(1,AllFrameStats.Orientation);
cA=cat(1,AllFrameStats.ConvexArea) * size_scale * size_scale;
fA=cat(1,AllFrameStats.FilledArea) * size_scale * size_scale;
eD=cat(1,AllFrameStats.EquivDiameter) * size_scale;
s=cat(1,AllFrameStats.Solidity);
ex=cat(1,AllFrameStats.Extent);
P=cat(1,AllFrameStats.Perimeter) * size_scale;
Pold=cat(1,AllFrameStats.PerimeterOld) * size_scale;
AllFrameConvPerims = cat(1,AllFrameStats.ConvexPerim) * size_scale;
AllFrameFL = cat(1,AllFrameStats.FiberLength) * size_scale;
AllFrameRadii = cat(1,AllFrameStats.MICRadii) * size_scale;
AllFrameBM = cat(1,AllFrameStats.BlebLength) * size_scale;
MICCenters = cat(1,AllFrameStats.MICCenter);
cx = MICCenters(:,1);
cy = MICCenters(:,2);
AllFrames = cat(1,AllFrameStats.Frame);


T=table(AllFrames,AllFrameNumbers, A, posx, posy, mjr, mnr, ecc, ori, ...
    cA, fA, eD, s, ex, P, Pold, AllFrameConvPerims, AllFrameFL, ...
    AllFrameRadii, AllFrameBM, cx, cy);

T.Properties.VariableNames={'Frames','CellNum','Area','CentroidX',...
    'CentroidY','MjrAxisLength','MnrAxisLength','Eccentricity',...
    'Orientation','ConvexArea','FilledArea','EquivDiameter',...
    'Solidity','Extent','Perimeter','PerimOld', 'ConvexPerim',...
    'FibLen', 'InscribeR','BlebLen','CenterX','CenterY'};

T=sortrows(T,2,'ascend');
idcs=strfind(save_path,'/');
writetable(T,[save_path(1:idcs(end)),save_path(idcs(end-1)+1:idcs(end)-1),'.csv'],'Delimiter',',','WriteRowNames',true)


%save the data just in case
save(sprintf([save_path,filesep,'workspace_checkpoint_%d.mat'],im_num),'AllFrameStats','AllFrameNumbers','blackoutlist')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm

disp('Complete :)')
disp('Please run Detect_Errors.m and ChangeCellNum.m if necessary')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm

disp('Complete :)')
disp('Please run Detect_Errors.m and ChangeCellNum.m if necessary')

function [Centers,Rs] = Max_Fit_Circles(BW)
    objs = unique(BW);
    Centers = zeros(length(objs)-1,2);
    Rs = zeros(length(objs)-1,1);
    for bb=2:length(objs)
        B=bwboundaries(BW==objs(bb),8,'noholes');
        grain=false(size(BW));
        for uu=1:length(B{1})
            grain(B{1}(uu,1),B{1}(uu,2))=true;
        end
        grain=uint8(grain)*255;
        [R,cx,cy] = max_inscribed_circle(grain,[]);
        Centers(bb-1,1)=cx;
        Centers(bb-1,2)=cy;
        Rs(bb-1) = R;
        if R<=3
            disp('SMALL RADIUS pause')
        end
    end
end



function ConvP = Convex_Perimeter(BW)
    %the bw image should already be filtered by size.
    %from black and white, create concave objects.
    bw_labeled=BW;
    ConvP=zeros(length(unique(bw_labeled))-1,1);
    %we might think that doing bwconvhull(bw,'objects') is indeed faster,
    %but it fails when two objects are nearby and when turned convex they
    %overlap.
    dnums = unique(bw_labeled);
    for i=2:length(dnums)
        objs = bwconvhull(bw_labeled==dnums(i));
        %calculate the perimeter.
        CC = bwconncomp(objs>0,8);
        stats = regionprops(CC, {'Perimeter'});
        ConvP(i-1)=stats.Perimeter;
    end
end

function FL = Fiber_Length(BW)
    %the bw image should already be filtered by size and in correct binary
    %form.
    bw2=BW;
    objs = unique(bw2);
    FL = zeros(length(objs)-1,1);
    for obj=2:length(objs)
        bw1 = bw2==objs(obj);
        %calculate longest skeleton path, aka fiber length for each object
        bw_skel=bwmorph(bw1,'skel','Inf');
        %how many objects are there?
        ep = find(bwmorph(bw_skel,'endpoints'));
        fl=0;
        for d=1:length(ep)
            D=bwdistgeodesic(bw_skel,ep(d),'quasi');
            if max(max(D))>fl && ~isnan(max(max(D)))
                fl=max(max(D));
            end
        end
        FL(obj-1) = fl;
        if isinf(fl)
            disp('INF FL pause')
        end
        if fl==0
            disp('ZERO FL pause')
        end
    end
end

%function newcellnums=FindPreviousCell3D(oldpos,oldcellnums,oldstats,oldz,pos,cellnums,stats,z)
%    %use change in shape 

function max_bleb_len = Bleb_Length(BW,Centers,Rs)
    bwl = BW;
    objs = unique(bwl);
    max_bleb_len = zeros(length(objs)-1,1);
    for obj=2:length(objs)
        %find skeleton of obj
        bw = bwl==objs(obj);
        bwskel = bwmorph(bw,'skel','Inf');
        Itemp = zeros(size(BW));
        Itemp(Centers(obj-1,2),Centers(obj-1,1))=1;
        Id = bwdist(Itemp);
        bw2 = bwskel;
        bw2(Id<ceil(Rs(obj-1)))=0;
        %find endpoints
        ep = find(bwmorph(bw2,'endpoints'));
        max_d = 0;
        for endpoint = 1:length(ep)
            bw3 = bwdistgeodesic(bw2,ep(endpoint),'quasi');
            m = max(bw3(~isinf(bw3)),[],'omitnan');
            if m > max_d
                max_d = m; 
            end
        end
        max_bleb_len(obj-1) = max_d;
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
