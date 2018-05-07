%% q1
clearvars
close all

% x and y axises limit from 0 to x_max and 0 to y_max respectively.
x_max = 100; %;
y_max = 100; %;

% distance of moving at every step
EPS = 2;

% maximum iterations
numNodes = 5000; %;

% attributions of starting point
q_start.coord = [0 0];
q_start.cost = 0;
q_start.parent = 0; %.parent means the index of parent node
 
% initialize the tree
nodes(1) = q_start;

% plot the goal area
figure(1)
axis([0 x_max 0 y_max])
goal_area = rectangle('Position',[70,45,5,5],'FaceColor',[0 .5 .5]);
xlabel('x')
ylabel('y')
hold on

%% grow the tree

for i = 1:1:numNodes
    
    % generate the random points in the given safe area and plot the points
    q_rand = [floor(rand(1)*x_max) floor(rand(1)*y_max)];
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    
    % Find the nearest point existing on the tree to the random point
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = norm(n.coord - q_rand);
        ndist = [ndist tmp];
    end
    [mini_distance, idx] = min(ndist);
    q_nearest = nodes(idx);
    
    % move to the random point with distance of eps if distance between
    % random point and nearest point is bigger than eps.
    q_new.coord = steer(q_rand, q_nearest.coord, mini_distance, EPS);
    line([q_nearest.coord(1), q_new.coord(1)], [q_nearest.coord(2), q_new.coord(2)],...
        'Color', 'k', 'LineWidth', 2);
    drawnow
    hold on
    q_new.cost = norm(q_new.coord - q_nearest.coord) + q_nearest.cost;
    q_new.parent = idx;
    
    % Append to nodes
    nodes = [nodes q_new];
    
    % Break if the link from second to last node to last node intersects any of
    % the four edges of the goal area
    if ~noCollision(q_nearest.coord, q_new.coord, [70,45,5,5])
        break
    end
end
q_end = q_new;
num_node_path = 1;
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)],...
        'Color', 'r', 'LineWidth', 2);
    hold on
    q_end = nodes(start);
    num_node_path = num_node_path+1;
end

%% total number of node in the tree
num_node_tree = length(nodes)

%% number of nodes in the sequence that reaches goal area
num_node_path

