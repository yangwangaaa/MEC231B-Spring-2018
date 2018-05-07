%% Q2
clearvars
close all

% x and y axises limit from 0 to x_max and 0 to y_max respectively.
x_max = 1000;
y_max = 7;

% distance of moving at every step
EPS = 2;

% maximum iterations
numNodes = 6000;

% attributions of starting point
q_start.coord = [0 2];
q_start.cost = 0;
q_start.parent = 0; %.parent means the index of parent node

% initialize the tree
nodes(1) = q_start;

% plot the safe area
figure(1)
x=[0 300 300 500 500 1000 1000 600 600 200 200 0]; %x coordinates of all the vertices
y=[0 0 4 4 0 0 3 3 7 7 3 3];  %y coordinates of all the vertices
X=[x,x(1)];   %????????????????????
Y=[y,y(1)];   %??
plot(X,Y,'k')  %?????
fill(x,y,'r')  % fill the safe zone with color
hold on

% plot the goal area
X_goal = [900,900,950,950,900];
Y_goal = [1,1.5,1.5,1,1];
plot(X_goal,Y_goal,'y')
fill(X_goal,Y_goal,'y')
%% grow the tree

for i = 1:1:numNodes
    pan = 0;
    % generate the random points in the given safe area and plot the points
    while ~pan
        q_rand = [rand*x_max,rand*y_max];
        pan = inpolygon(q_rand(1),q_rand(2),X,Y);
    end
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = norm(n.coord - q_rand);
        ndist = [ndist tmp];
    end
    [mini_distance, idx] = min(ndist);
    q_nearest = nodes(idx);
    q_new.coord = steer(q_rand, q_nearest.coord, mini_distance, EPS);
    
    InorOn = inpolygon(q_new.coord(1),q_new.coord(2),X,Y);
    if InorOn
        line([q_nearest.coord(1), q_new.coord(1)], [q_nearest.coord(2), q_new.coord(2)],...
            'Color', 'k', 'LineWidth', 2);
        drawnow
        hold on
        q_new.cost = norm(q_new.coord - q_nearest.coord) + q_nearest.cost;
        q_new.parent = idx;
        
        % Break if the link from second to last node to last node intersects any of
        % the four edges of the goal area
        nodes = [nodes q_new];
    end
    
    % move to the random point with distance of eps if distance between
    % random point and nearest point is bigger than eps.
    
    if ~noCollision(q_nearest.coord, q_new.coord, [900,1,50,0.5])
        break
    end
end

q_end = q_new;
num_node_path = 1;
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)],...
        'Color', 'g', 'LineWidth', 2);
    hold on
    q_end = nodes(start);
    num_node_path = num_node_path+1;
end

%% total number of node in the tree
num_node_tree = length(nodes)

%% number of nodes in the sequence that reaches goal area
num_node_path