
%% THIS IS FOR OBJECT II TO MEO
%% Define inputs

% Crazy state at burn 1
state1 = [1.096035657606564e+04;2.698333735244205e+04;1.194164375722560e+02;-2.985109599699163;-1.186895083866727;0.014082710185299];

% MEO state at arrival (burn 2)
state2 = [-1.783230987211009e+04;1.214512502799535e+04;-4.006338669676189e+02;-3.772645841393155;-3.030358641398373;-0.042496438145642];

deltat = 7.113944104937772e+03; % seconds

v1 = [-4.13766045630268; -0.0169025030537010; -0.0825223428910495]; % from Shaniya
v2 = [-2.83802067344737; -4.31769720933036; -0.0407483702417452];
% MAKE SURE V2 IS THE SAME AFTER ODE


%% Prepare to call ode

options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
tstart = 0;
tstop = deltat; % sec
tspan = [tstart tstop];

r = state1(1:3);
v = v1(1:3);
state = [r;v];

% Call function that does the work
[newtime,newstate] = ode45(@non_impulsive_COAST,tspan,state,options,mu);

%% Plot it

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % S/C at Burn 1
       plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2) % orbit traj
       plot3(state1(1),state1(2),state1(3),'*','LineWidth',2) 

   % S/C at Burn 2
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) % orbit traj
       plot3(state2(1),state2(2),state2(3),'*','LineWidth',2) % current location

   % Lamberts Traj
       plot3(newstate(:,1),newstate(:,2),newstate(:,3),'k','LineWidth',2) % orbit traj

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Phase III Lamberts Transfer']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 9) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','Object II Orbit','S/C at Burn 1', 'MEO1 Orbit','S/C at Burn 2','Lamberts Trajectory','location','best')



        %% THIS IS FOR MEO TO GEO
%% Define inputs

% MEO state at burn 1
state1 = [-1.473446403861598e+04;-3.762629130736880e+04;28.673587265657297;2.189260581538516;-1.186591955263678;0.046887766408053];

% GEO state at arrival (burn 2)
state2 = [3.433878743797304e+04;-2.417490245087891e+04;-1.033902024430127e+03;1.771289141948836;2.524834624894050;-0.048552577088648];

deltat = 1.258729719024421e+04; % seconds

v1 = [4.01987248995049; -0.531845543271505; -0.101165539178234];
v2 = [3.12988713997306; 2.42945783150817; -0.0541848141634415];

% MAKE SURE V2 IS THE SAME AFTER ODE


%% Prepare to call ode

options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
tstart = 0;
tstop = deltat; % sec
tspan = [tstart tstop];

r = state1(1:3);
v = v1(1:3);
state = [r;v];

% Call function that does the work
[newtime,newstate] = ode45(@non_impulsive_COAST,tspan,state,options,mu);

%% Plot it

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % S/C at Burn 1
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) % orbit traj
       plot3(state1(1),state1(2),state1(3),'*','LineWidth',2) 

   % S/C at Burn 2
       plot3(geo1.state_t0(:,1),geo1.state_t0(:,2),geo1.state_t0(:,3),'k','LineWidth',2) % orbit traj
       plot3(state2(1),state2(2),state2(3),'*','LineWidth',2) % current location

   % Lamberts Traj
       plot3(newstate(:,1),newstate(:,2),newstate(:,3),'r','LineWidth',2) % orbit traj

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Phase IV Lamberts Transfer']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 9) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','MEO1 Orbit','S/C at Burn 1', 'GEO1 Orbit','S/C at Burn 2','Lamberts Trajectory','location','best')
