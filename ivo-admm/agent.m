classdef agent 

    properties
        index               % agent index
        radius 
        pose
        vel
        vel_max

        goal

        r_obstacles
        pose_obstacles
        vel_obstacles
        
        r_neighbours
        pose_neighbours
        vel_neighbours
        
        safety_margin       % collision safety margin
        
        n                   % number of active neighbours for agent

        u                   % control

        z                   % current estimate of control
        x                   % estimate of contol of neighbours
        y                   % dual variable associated with constraint of each neighbouring agent 
        xN                  % estimate of agents control that neighbour has
        yN                  % dual varibale communicated by neighbour
        zN                  % communicated neighbours control
        
        rho                 % multiplier for augmented square term     
        alpha               % smoothness parameter
        beta                % paramter for agent-agent constraints
        lambda_range        % range of dual variable when solving z-update
        avoidance_range     % obstacle avoidance range
    end

    methods 
        function default_agent = agent(env, index)
            % creats a default agent based on the 
            default_agent.index = index;
            default_agent.radius = env.agents_radius;
            default_agent.pose = env.agents_pose(index,:);
            default_agent.vel = [0 0];
            default_agent.vel_max = 2;

            default_agent.goal = env.agents_goal(index,:);
            
            default_agent.r_obstacles = env.obstacles_radius; % radius of all obstacles that are in the env
            default_agent.pose_obstacles = env.obstacles_pose;
            default_agent.vel_obstacles = env.obstacels_vel;

            default_agent.r_neighbours = repmat(env.agents_radius, env.n_agents, 1); % radius of all agents that are in the env
            pose_neighbours = env.agents_pose;                % pose of other agents acting as obstacles
            default_agent.pose_neighbours = pose_neighbours; 

            vel_neighbours = env.agents_vel;
            default_agent.vel_neighbours = vel_neighbours;    % pose of other agents acting as obstacles

            default_agent.safety_margin = 2*default_agent.radius;
            default_agent.n = env.n_agents-1;               % considering all neighbours to be active

            default_agent.u = [0 0]; % (default_agent.goal-default_agent.pose)*eps;

            default_agent.z = [0 0]; % (default_agent.goal-default_agent.pose)*eps;
            default_agent.x = zeros(env.n_agents,2);
            default_agent.y = zeros(env.n_agents,2);

            default_agent.xN = zeros(env.n_agents,2);
            default_agent.yN = zeros(env.n_agents,2);
            default_agent.zN = env.agents_zN;

            for i=1:env.n_agents    % retrieve x and y from the environment (ideally other agents)
                agent.xN(i,:) = env.agents_xN(i,:,i); 
                agent.yN(i,:) = env.agents_yN(i,:,i);
            end
            
            default_agent.rho = 0.5;
            default_agent.alpha = 1;
            default_agent.beta =  0.365./vecnorm(env.agents_pose - default_agent.pose,2,2).^2; default_agent.beta(index) = 0;
            % default_agent.beta = 0*default_agent.beta; default_agent.beta(index) = 0; % testing with beta = 0
            % default_agent.beta = 1*[1; 1]; default_agent.beta(index) = 0;
            default_agent.lambda_range = [0 1e10];

            default_agent.avoidance_range = 6;
        end

        function [new_agent, new_env] = z_update(agent,env)

            for i=1:env.n_agents    % retrieve x and y from the environment (ideally other agents)
                agent.xN(i,:) = env.agents_xN(i,:,i); 
                agent.yN(i,:) = env.agents_yN(i,:,i);
            end
            
            vel_desired = agent.vel_max*(agent.goal-agent.pose)/norm(agent.goal-agent.pose);
            e = -( 2*(vel_desired-agent.vel) +  sum(agent.yN + agent.rho*agent.xN) )' / (1 + agent.alpha + 0.5*agent.rho*(agent.n+1));
            
            p_relative = (agent.pose_obstacles-agent.pose)'; v_relative = (agent.vel_obstacles-agent.vel)'; 
            P = p_relative*p_relative' + (agent.r_obstacles+agent.radius+agent.safety_margin)^2 * eye(2);
            q = -2*( P*v_relative + (norm(p_relative)^2) * (agent.z'-v_relative) );
            r = v_relative'*P*v_relative + (norm(p_relative)^2) * (norm(agent.z'-v_relative)^2) + 2*norm(p_relative)^2 * (agent.z'-v_relative)'*v_relative;

            % neglecting obstacle
            % P = 0*P; q = 0*q; r = 0;

            % formulate the problem for z-update
            problem.e = e;
            problem.P = P;
            problem.q = q;
            problem.r = r;
            problem.z0 = agent.z';
            problem.lambda_range = agent.lambda_range; problem.index = agent.index;

            z_hat = fminConSolveQCQP(problem); 
            % z_hat = cvxSolveQCQP(problem); 
            % z_hat = solveQCQP(problem);

            new_agent = agent; 
            new_agent.z = z_hat';

            new_env = env; 
            new_env.agents_zN(agent.index,:) = z_hat'; 
        end
    
        function [new_agent, new_env] = x_update(agent,env)
            
            idx = agent.index;
            agent.zN = env.agents_zN; % retrieve neighbours current control estimate

            A = zeros(2*env.n_agents); b = zeros(2*env.n_agents,1);
            p_relative = agent.pose_neighbours - agent.pose;
            v_relative = agent.vel_neighbours - agent.vel;
            for j = 1:agent.n+1             % iterating over all neighbours and self
                if j ~= idx
                    Dj = zeros(2,2*(agent.n+1)); Dj(:,2*j-1:2*j) = eye(2); Dj(:,2*idx-1:2*idx) = -eye(2);
                else
                    Dj = zeros(2,2*(agent.n+1));
                end
                Cj = zeros(2,2*(agent.n+1)); Cj(:,2*j-1:2*j) = eye(2);  
                pj = p_relative(j,:)'; vj = v_relative(j,:)';
                A = A + agent.beta(j)*Dj'*(pj*pj' + (agent.r_neighbours(j) + agent.r_neighbours(idx) + agent.safety_margin)^2 * eye(2))*Dj + 0.5*agent.rho*(Cj'*Cj);
                b = b + agent.beta(j)*Dj'*(( (pj*pj' + ((agent.r_neighbours(j) + agent.r_neighbours(idx) + agent.safety_margin)^2 - norm(pj)^2) * eye(2)) )*vj - norm(pj)^2 * Dj * reshape(agent.zN',2*env.n_agents,1)) + 0.5*Cj'*(agent.y(j,:)-agent.rho*agent.zN(j,:))';
            end
            new_agent = agent;
            new_agent.x = reshape(-A\b,2,env.n_agents)'; new_agent.xN(idx,:) = new_agent.x(idx,:); % solve for x

            fprintf('x = [%f %f]\n',new_agent.x(idx,1),new_agent.x(idx,2));

            new_env = env; new_env.agents_xN(:,:,idx) = new_agent.x; % communicate x to the environment
        end

        function [new_agent, new_env] = y_update(agent,env)
            new_agent = agent;
            new_agent.y = agent.y + agent.rho*(agent.x - agent.zN); % all y for an agent will be same if they start from same initial value

            new_env = env; new_env.agents_yN(:,:,agent.index) = new_agent.y;
        end

        function [c,ceq] = generate_constraints(agent,u)
            % generate constraints
            n_obs = size(agent.pose_neighbours,1);
            c = zeros(1,n_obs);
            for i = 1:n_obs
                p_relative = (agent.pose_neighbours(i,:)-agent.pose); v_relative = (agent.vel_neighbours(i,:)-agent.vel); 
                if p_relative*v_relative' <= norm(p_relative)*norm(v_relative)*cosd(80) && norm(p_relative) <= agent.avoidance_range
                    critical_dist = norm(p_relative)^2 - (agent.radius + agent.safety_margin + agent.r_obstacles(i))^2;
                    c(i) = ((p_relative*(v_relative-u)')/norm(v_relative-u))^2 - critical_dist;
                else
                    c(i) = 0;
                end
                ceq = [];
            end
        end

        function [new_env, new_agent] = move_and_update(agent,env)
            
            p_obstacles = env.agents_pose;
            p_obstacles(agent.index,:) = [];
            agent.pose_neighbours(1:env.n_agents-1,:) = p_obstacles;

            v_obstacles = env.agents_vel;
            v_obstacles(agent.index,:) = [];
            agent.vel_neighbours(1:env.n_agents-1,:) = v_obstacles;
            
            if norm(agent.pose-agent.goal) > 1e-1
                vel_desired = agent.vel_max*(agent.goal-agent.pose)/norm(agent.goal-agent.pose);
                objective = @(u) (norm(vel_desired-(agent.vel+u)))^2 + agent.alpha * (norm(u)^2);
                constraints = @(u) agent.generate_constraints(u);

                problem.objective = objective;
                problem.nonlcon = constraints;

                options = optimoptions('fmincon','display','off','Algorithm','sqp');
                problem.options = options;
                problem.solver = 'fmincon';
                problem.x0 = agent.u;
                
                agent.u = fmincon(problem);
                
                v = agent.vel + agent.u; 
                if norm(v)>agent.vel_max
                    agent.vel = agent.vel_max*v/norm(v);
                else
                    agent.vel = v;
                end
            else
                agent.vel = 0;
                agent.u = 0;
            end

            agent.pose = agent.pose + agent.vel*env.deltaT;

            new_env = env;
            new_env.agents_pose(agent.index,:) = agent.pose;
            new_env.agents_vel(agent.index,:) = agent.vel;

            new_agent = agent;
        end

        function [new_env, new_agent] = move(agent,env)
            if norm(agent.pose-agent.goal) > 1e-1
                agent.u = agent.z;
                v = agent.vel + agent.u;
                if norm(v)>agent.vel_max
                    agent.vel = agent.vel_max*v/norm(v);
                else
                    agent.vel = v;
                end
            else
                agent.vel = 0;
                agent.u = 0;
            end

            agent.pose = agent.pose + agent.vel*env.deltaT;

            new_env = env;
            new_env.agents_pose(agent.index,:) = agent.pose;
            new_env.agents_vel(agent.index,:) = agent.vel;

            new_agent = agent;
        end

        function new_agent = reinit(agent)
            % reinitialize certain variables
            new_agent = agent;
            new_agent.z = [0 0];
            
            new_agent.x = 0*new_agent.x;
            new_agent.y = 0*new_agent.y;

            new_agent.xN = 0*new_agent.xN;
            new_agent.yN = 0*new_agent.yN;
            new_agent.zN = 0*new_agent.zN;
        end

    end
end


% cvx-solver for z-update
function z = cvxSolveQCQP(problem)
    P = problem.P; q = problem.q; r = problem.r; e = problem.e;
    cvx_begin quiet
        variable z(2)
        dual variables lam
        minimize( z'*z + e'*z )
        lam: z'*P*z + q'*z + r <= 0;
    cvx_end
end

% fmincon-solver for z-update
function z = fminConSolveQCQP(problem)
    P = problem.P; q = problem.q; r = problem.r; e = problem.e; z0 = problem.z0;
    fun = @(z)quadobj(z,e);
    quadconstr = @(z)quadcon(z,P,q,r);
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',@(z,lambda)quadhess(z,lambda,P),'Display','off');
    try
        [z,fval,eflag,output,lambda] = fmincon(fun,z0,[],[],[],[],[],[],quadconstr,options);
    catch
        keyboard
    end
    fprintf('Agent %d:\n',problem.index)
    fprintf('lambda %f \nz = [%f %f] \n',lambda.ineqnonlin,z(1),z(2));
end

function [y,grady] = quadobj(x,e)
    y = x'*x + e'*x;
    if nargout > 1
        grady = 2*x + e;
    end
end

function [y,yeq,grady,gradyeq] = quadcon(x,P,q,r)
    y = x'*P*x+q'*x+r; yeq = [];
    grady = 2*P*x + q; gradyeq = [];
end

function hess = quadhess(x,lambda,P)
    hess = 2*(eye(2)+lambda.ineqnonlin*P);
end

% solver for z-update
function z = solveQCQP(problem)
    P = problem.P; q = problem.q; r = problem.r; e = problem.e;
    lambda_range = problem.lambda_range; z0 = problem.z0;
    lambda = lambda_range(1);

    while(diff(lambda_range)>1e-3)
        z_hat = -0.5*(eye(2)+lambda*P)\(e+lambda*q);

        if z_hat'*P*z_hat + q'*z_hat + r > 0
            lambda_range(1) = lambda;
        else
            lambda_range(2) = lambda;
        end
        lambda = mean(lambda_range);
    end

    fprintf('lambda: %f\n',lambda);

    if z_hat'*P*z_hat + q'*z_hat + r <= 0
        fprintf('Feasible solution found for z-update\n');
        z = z_hat;
    else
        fprintf('Solution infeasible for z-update\n');
        z = z_hat;
    end
end
