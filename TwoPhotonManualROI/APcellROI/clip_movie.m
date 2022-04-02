function  clip_movie(fn, fr_range)
% load fn tiff movive
% replace it by subselecting the range [start stop_frames]
% use Aki's functions

%INPUT:
%   fn: fullfpath of the movie
%   fr_range: [start stop]
%%

[stack,info,~] = read_tiff(fn);
write_tiff(fn,stack(:,:,[fr_range(1):fr_range(2)]),info);