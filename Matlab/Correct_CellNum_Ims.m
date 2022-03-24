


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We want to ask the user to identify the max projected stack of images
%which we will use to identify the cells. 
disp('Select max-projected stack of images  ')
gray_path=uigetdir(pwd,'select max-projected stack of images');
gray_ims=dir([gray_path,filesep,'*png']);
gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files

if length(gray_ims)==0
    gray_ims=dir([gray_path,filesep,'*tif']);
    gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select where you want to save the printed images to.
disp('Select directory you wish to save printed images to ')
save_path=uigetdir(pwd,'Select directory you wish to save printed images to');
if strcmp(save_path,gray_path)
    error('You do not want gray path and save path to be the same')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




disp("Please select .CSV file that was saved after running cell tracking: ")

[file,path]=uigetfile('*.csv');

T=readtable(fullfile(path,file));


%F = %unique(T.Frames);

for frame = 1:length(gray_ims)
    OneFrameCellnums = T.CellNum(T.Frames==frame);
    xs = T.CenterX(T.Frames==frame);
    ys = T.CenterY(T.Frames==frame);
    OneFrameCenters=[xs,ys];
    Igray=im2uint8(imread([gray_path,filesep,gray_ims(frame).name]));
    figure
    %Igray=imread([gray_path,filesep,gray_ims(im_num).name]);
    imshow(Igray)
    if ~isempty(OneFrameCellnums)
        for k=1:size(OneFrameCenters,1)
            %if sum(blackoutlist==OneFrameCellnums(k))==0
            c=OneFrameCenters(k,:);
            text(c(1),c(2),sprintf('%d',OneFrameCellnums(k)),...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle',...
                'FontWeight','bold',...
                'Color','r');
            %end
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
    imwrite(im,[save_path,filesep,gray_ims(frame).name])
end