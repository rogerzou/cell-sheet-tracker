
format short g;
format compact;

%% load images

GT = imreadgroundtruth('img/DDC1_GT.png', true);  % loads ground truth of first frame
ALL = loadtiff('img/DDC1_all.tif');               % loads tiff stack into 3D image matrix
xrange = 50:150; yrange = 50:150; zrange = 1:10;  % crop data set
ALL = ALL(yrange, xrange, zrange);
GT = GT(yrange, xrange);
GT = padarray(GT, [20,20]);                       % zero padding to avoid crashing due to sampling outside image
ALL = padarray(ALL, [20,20,0]);

%% obtain processed images (smoothed and segmentations)

N = size(ALL,3);                    % number of images

% gaussian smoothed images in ALL_s
H = fspecial('gaussian', 5, 2);
ALL_s = ALL;                                                                % ALL_s holds the smoothed image stack
for ii=1:N
    ALL_s(:,:,ii) = imfilter(ALL(:,:,ii), H);                               % filter each image frame
end

% watershed segmentations of all image frames in SEG
hminima = 25;                                                               % the h parameter for matlab's imhmin function optimized to this data set
SEG = false(size(ALL));                                                     % SEG holds the segmentations of image stack
for ii=1:N
    L = watershed( imhmin(medfilt2(ALL(:,:,ii),[3,3]), hminima) );          % watershed segmentation on each frame (after light median filtering for noise removal)
    SEG(:,:,ii) = imreadgroundtruth(L==0, true);                            % convert segmentation to ground truth format
end

% NOTE: to use a segmented first frame rather than the default ground truth, uncomment the following line:
% GT = SEG(:,:,1);

%% optimization

% algorithm parameters
l = 17;                 % width of edge window
w = 25;                 % width of vertex window
alpha = 0.5;            % scaling edge cost contributions
interval = 1;           % merge vertices if distance is less or equal
spacing = 15;           % spacing of control points on splines
parallel = false;       % parallelization with parfor
verboseE = 0;           % verbose flag for edge optimization
verboseG = 0;           % verbose flag for vertex optimization
siftflow = true;        % SIFT flow flag
fname = ['tmp_grid2_', datestr(clock, 'mmddyy_HH:MM:SS')];
options = struct('l',l,'w',w,'alpha',alpha,'interval',interval, ...
    'spacing',spacing,'parallel',parallel,'verboseE',verboseE, ...
    'verboseG',verboseG,'siftflow',siftflow,'fname',fname);

% track graph
t1 = tic;
data = membraneTrack(ALL, GT, options);     % this is the function that performs the tracking
disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
fprintf('\nDone!\n');
toc(t1);

%% clean up

% save relevant variables
save('z_data.mat', 'data');  % data is a struct that stores all raw tracing information

% create tracking movie
Movie = [];
for ii=1:size(ALL,3)
    fig = displayGraph(ALL(:,:,ii), data(ii).VALL, data(ii).EALL, 'on');    % function parses the data struct to display tracks over original image
    Movie = [Movie, immovie(print(fig, '-RGBImage'));];                     % save current frame to movie
    close(fig);
end
fnameall = 'z_movie.avi';
writerObj = VideoWriter(fnameall);  % create VideoWriter object to write frames to AVI movie
writerObj.FrameRate = 2; writerObj.Quality = 100;
open(writerObj); writeVideo(writerObj, Movie); close(writerObj);
