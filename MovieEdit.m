%% 19 Jan, 2021
% More details in
% https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html
%% Create the Image Sequence
workingDir = 'G:\Judy\';
mkdir(workingDir)
mkdir(workingDir, 'images')

shuttleVideo = VideoReader('G:\Judy\ViewAmp 2021-01-19 10-19-30.mp4');

ii = 1;

while hasFrame(shuttleVideo)
   img = readFrame(shuttleVideo);
   filename = [sprintf('%05d', ii) '.jpg'];
   fullname = fullfile(workingDir, 'images', filename);
   imwrite(img, fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   ii = ii + 1;
end

%% Create New Video with the Image Sequence
imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'ViewAmp.avi')); % Prepare a avi', 'mp4', or 'mj2' file
outputVideo.FrameRate = 24; % The same as the fps
open(outputVideo);
for ii = 1:72:length(imageNames)  % The frame number for video
    frame = imread(fullfile(workingDir,'images',imageNames{ii}));
    writeVideo(outputVideo, frame);
end
close(outputVideo);