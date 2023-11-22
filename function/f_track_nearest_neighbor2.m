function [ matches ] = f_track_nearest_neighbor2( part_A, part_B, f_o_s )
%FUNCTION matches = f_track_nearest_neighbor2(part_A, part_B, f_o_s)
% Objective: Nearest neigbhor tracking routine.
% ---------------------------------------------------
%
%   INPUT:      part_A      - coordinates of particles in image A [n x 2]
%    
%               part_B      - coordinates of particles in image B [n x 2]
% 
%               f_o_s       - field of search [px]
%
%   OUTPUT:     matches     - list of indices of matching particles [m x 2]
%    
% ---------------------------------------------------
% References
% [1] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% ---------------------------------------------------
% origin: Thomas Janke / 2017.09
% modified: Jin Yang / 2020.12, Alex Landauer / 2023.11
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
for parInd = 1:size(part_A,1) % Loop all particles in Image A
      
    temp_dist = sqrt((part_B(:,1)-part_A(parInd,1)).^2 + (part_B(:,2)-part_A(parInd,2)).^2  ); 
    % Calculate all possible distances between particle parInd and particles in image B
    
    if min(temp_dist) < f_o_s % Just consider candidates within field of search
        [~,index_min] = min(temp_dist); % Find candidate with minimum displacement
        matches{parInd} = [parInd, index_min];
    end
         
end

% Check if there are found matches
if exist('matches','var')
    matches = cell2mat(matches');
else
    matches = [];
end

end

