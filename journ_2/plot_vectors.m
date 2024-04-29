function plot_vectors(vectors,normalize)
% 
% PLOT_VECTORS plots vectors starting at the origin. Useful to compare sets
% of vectors.
% 
% 
%USAGE
%-----
% plot_vectors(vectors)
% plot_vectors(vectors,normalize)
% 
% 
%INPUT
%-----
% - VECTORS  : cell array with Vix2 or Vix3 matrices (i=1,2,...) in the
%   cell i, with the x-, y- and z-coordinates for Vi vectors
% - NORMALIZE: 1 (normalize the vectors) or 0 (don't normalize)
% 
% 
%OUTPUT
%------
% Figure with the vectors from VECTORS{i} starting at the origin. Colors
% for VECTORS{i}: blue (i=1,8,...), red (i=2,9,...), green (i=3,10,...),
% cyan (i=4,11,...), magenta (i=5,12,...), yellow (i=6,13,...), black
% (i=7,14,...).
% 
% 
% See also COMPASS, QUIVER, QUIVER3, PLOT3

% Guilherme Coco Beltramini (guicoco@gmail.com)
% 2012-Jul-20, 10:21am


vec_color = ['b' 'r' 'g' 'c' 'm' 'y' 'k'];
num_color = length(vec_color);


% Input
%==========================================================================
if nargin<2
    normalize = 0;
end
if ~iscell(vectors)
    vectors = {vectors};
end
if ~isequal(normalize,0) && ~isequal(normalize,1)
    error('Unknown option for "normalize"')
end

vectors1=readmatrix('s.dat');
% Plot figure
%==========================================================================
figure
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p, 'FaceColor', 'r')
axis equal
hold on

for vv=1:length(vectors)
    
    dim = size(vectors{vv},2);
    if dim~=2 && dim~=3
        fprintf('Vector %d ignored: the vectors must be Vx2 or Vx3 matrices.\n',vv)

        continue
    end
    
    Vsz = size(vectors{vv},1);
    if normalize
        for gg=1:Vsz
            vectors{vv}(gg,:) = vectors{vv}(gg,:)/norm(vectors{vv}(gg,:)');
        end
    end
    if dim==2
        

        quiver(vectors1(:,1),vectors1(:,2),vectors{vv}(:,1),vectors{vv}(:,2),5,'Color',vec_color(rem(vv-1,num_color)+1),'LineWidth',1)
   

    elseif dim==3
        quiver3(zeros(Vsz,1),zeros(Vsz,1),zeros(Vsz,1),...
            vectors{vv}(:,1),vectors{vv}(:,2),vectors{vv}(:,3),...
            'Color',vec_color(rem(vv-1,num_color)+1),'LineWidth',1.5)
    end
    
end


% Axes
%--------------------------------------------------------------------------
if normalize
    plot3([-1 1],[ 0 0],[ 0 0],'-k','LineWidth',1)
    plot3([ 0 0],[-1 1],[ 0 0],'-k','LineWidth',1)
    plot3([ 0 0],[ 0 0],[-1 1],'-k','LineWidth',1)
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end


% Labels
%--------------------------------------------------------------------------
xlabel('Vx')
ylabel('Vy')
zlabel('Vz')

grid on
hold on