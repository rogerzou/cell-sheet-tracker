
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
