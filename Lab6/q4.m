%% Q4
clearvars
close all

% x and y axises limit from 0 to x_max and 0 to y_max respectively.
x_min = 0; x_max = 1000;
y_min = -20; y_max = 20;
psi_min = -pi; psi_max = pi; 

% distance of moving at every step
EPS = 2;

% maximum iterations
numNodes = 6000;

% attributions of starting point
q_start.coord = [0 0 0];
q_start.cost = 0;
q_start.parent = 0; %.parent means the index of parent node

% initialize the tree
nodes(1) = q_start;

% plot the safe area
figure(1)
x=[0 300 300 500 500 1000 1000 600 600 200 200 0]; %x coordinates of all the vertices
y=[0 0 4 4 0 0 3 3 7 7 3 3];  %y coordinates of all the vertices
X=[x,x(1)];   
Y=[y,y(1)];   
plot(X,Y,'k') 
fill(x,y,'r')  % fill the safe zone with color
hold on

% plot the goal area
X_goal = [900,900,950,950,900];
Y_goal = [1,1.5,1.5,1,1];
plot(X_goal,Y_goal,'y')
fill(X_goal,Y_goal,'y')
%% grow the tree
tic;
for i = 1:1:numNodes
    pan = 0;
    % generate the random points in the given safe area and plot the points
    while ~pan
        q_rand = [rand*(x_max-x_min)+x_min,rand*(y_max-y_min)+y_min,rand*(psi_max-psi_min)+psi_min];
        pan = inpolygon(q_rand(1),q_rand(2),X,Y);
    end
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = norm(n.coord(1:2) - q_rand(1:2));
        ndist = [ndist tmp];
    end
    [mini_distance, idx] = min(ndist);
    q_nearest = nodes(idx);
    
    % Get the new Point
    % q_new.coord = steer(q_rand, q_nearest.coord, mini_distance, EPS);
    % Instead of using steer function directly, we apply the brute force
    % approach.
    
    % Parameters
    vx = 30; L = 3;
    dxdt = @(t,x,delta) [vx*cos(x(3));vx*sin(x(3));vx/L*tan(delta)];
    direc_desired = (q_rand(1:2) - q_nearest.coord(1:2))./norm(q_rand(1:2) - q_nearest.coord(1:2));
    direc = 0;
    flag_fea = false;
    for delta = -20/180*pi:2/180*pi:20/180*pi
        safety = true;
        [t,x] = ode45(@(t,x) dxdt(t,x,delta),[0 0.1],q_nearest.coord');
        q_fea = x(end,:);
        direc_fea = (q_fea(1:2) - q_nearest.coord(1:2))./norm(q_fea(1:2) - q_nearest.coord(1:2));
        % to check if any step there is a collison
        stepNum = size(x,1);
        for i = 1:stepNum
            if ~inpolygon(x(i,1),x(i,2),X,Y)
                safety = false;
                break;
            end
        end
        
        if sum(direc_fea.*direc_desired) > direc && safety
            flag_fea = true;
            direc = sum(direc_fea.*direc_desired);
            q_new.coord = (q_fea - q_nearest.coord)./norm(q_fea - q_nearest.coord) * EPS + q_nearest.coord;
        end
    end
    InorOn = inpolygon(q_new.coord(1),q_new.coord(2),X,Y);
    if flag_fea && InorOn
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
    
    if ~noCollision(q_nearest.coord(1:2), q_new.coord(1:2), [900,-1,50,2]) 
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

%% Total Calculation Time
toc