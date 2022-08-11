clear
%tridiag(a,b,c,f)

% Values for Crank-Nicolson Grid Initialization
D = 4080;  %Diffusion Strength
tau = 7.5;
L = 2000/(sqrt(tau*D)); %Length of domain
M = 100; %Number of stripes/segments
M1 = 5000; %Number of spatial grid points
dx = L/M1; %Spatial Step
dt = .0001;  %Time Step
Time = 20; %Total Sim Time
Track = 0; %Track images
Delay = 10000; %Delay Time (Must be Integer indicating time step delay)
r = D*dt/(2.0*dx*dx); %CN Parameter

%Parameters
phase_shift_value = 0.25; % Corresponds to initial phase difference
na0p = 985.56;
na1p = 17991.82;
K1 = 5937;
n1 = 4;
K2 = 10;
n2 = 2;
de = 2257;
Ke = 1000; 
nr0p = 74.46;
nr1p = 46495.80;
ne = 4;
gamma = 0.128; %Natural Degradation Rate

% Non-Dimensionalized Parameters
na0p_s = na0p*tau/(K1); 
na1p_s = na1p*tau/(K1); 
nr0p_s = nr0p*tau/(K1); 
nr1p_s = nr1p*tau/(K1); 
D_s = 1;
K2_s = K2/K1;
de_s = de*tau;
Ke_s = Ke/K1;
gamma_s = gamma*tau;
x0 = sqrt(D*tau);
r_s = D_s*dt/(2.0*dx*dx);

% Initialize phases of subpopulations:
for i=1:M1
    q1(i,1) = 0; %Initial condition
    q2(i,1) = 0; %Initial condition
end

for i = 1:M1
    %Initial Delay Terms are zero
    q1_del(i) = 0.0;
    q2_del(i) = 0.0;
    %Setup Tridiagonal Matrix (Diffusion and Decay Terms)
    A(i) = -r_s;
    B(i) = 2.0*r_s+1.0+gamma_s*dt/2.0;
    C(i) = -r_s;
    if (i == 1) % Implement open boundary conditions (set q1,q2 = 0 at boundary)
        A(i) = 0.0;
        B(i) = 1.0;
        C(i) = 0.0;
        F1(i) = 0.0;
        F2(i) = 0.0;
    elseif (i == M1)
        A(i) = 0.0;
        B(i) = 1.0;
        C(i) = 0.0;
        F1(i) = 0.0;
        F2(i) = 0.0;
    end
end

%Piecewise define production terms so activator production is on for only
%even stripes and repressor production is on for only odd stripes.
for i = 1:M1
    if (mod(floor(i/(M1/M)),2) == 0) % CHANGE!
        na0_s(i) = na0p_s;
        na1_s(i) = na1p_s;
        nr0_s(i) = 0.0;
        nr1_s(i) = 0.0;
    else
        na0_s(i) = 0.0;
        na1_s(i) = 0.0;
        nr0_s(i) = nr0p_s;
        nr1_s(i) = nr1p_s;
    end
end

%Spatial Grid
for i = 1:M1
    x(i) = (i-1)*dx;
end

for k = 1:Time/dt
    for i = 2:M1-1
    %Compute RHS for both q1 and q2 each time step
    F1(i) = q1(i+1,k)*r_s+q1(i,k)*(1-2.0*r_s-dt*gamma_s/2.0)+q1(i-1,k)*r_s+ dt*(na0_s(i)+na1_s(i)*(q1_del(i))^(n1))/(1+(q1_del(i))^n1+(q2_del(i)/K2_s)^n2) - dt*de_s*(q1(i,k)*(q2_del(i)/Ke_s)^ne)/(1+(q2_del(i)/Ke_s)^ne);
    F2(i) = q2(i+1,k)*r_s+q2(i,k)*(1-2.0*r_s-dt*gamma_s/2.0)+q2(i-1,k)*r_s+ dt*(nr0_s(i)+nr1_s(i)*(q1_del(i))^(n1))/(1+(q1_del(i))^n1) - dt*de_s*(q2(i,k)*(q2_del(i)/Ke_s)^ne)/(1+(q2_del(i)/Ke_s)^ne);
    end
    
    % Solve the tridiagonal matrix to get the next time step values
    q1(:,k+1) = tridiag(B,A,C,F1);
    q2(:,k+1) = tridiag(B,A,C,F2);
  
    %Delay Term
    for i = 1:M1
        if (k <= Delay)
            q1_del(i) = 0.0;
            q2_del(i) = 0.0;
        else
            q1_del(i) = q1(i,k-Delay);
            q2_del(i) = q2(i,k-Delay);
        end
    end
    

    % PRINTS THE TIME STEP!
    if mod(k,500) == 0
        k
    end
    
    
end

% Record time and concentration values for midpoint.
sim_time = 1:Time/dt+1;
midpoint_activator_concentration = q1(M1/2,:);
midpoint_repressor_concentration = q2(M1/2,:);

% Get normalized concentrations
norm_midpoint_activator_concentration = midpoint_activator_concentration/max(midpoint_activator_concentration);
norm_midpoint_repressor_concentration = midpoint_repressor_concentration/max(midpoint_repressor_concentration);

% Get minima for midpoint q1 concentrations to extract phase:
TF = islocalmin(q1(M1/2,:));
q1(M1/2,TF);
min_times = sim_time(TF);

% Plot normalized midpoint activator and repressor concentration vs. time
% along with minima so you can evaluate:
clf()
fig1 = figure(1);
clf()
set(fig1,'visible','on');
hold on
gca()
plot(sim_time,norm_midpoint_activator_concentration,'ob-')
plot(sim_time,norm_midpoint_repressor_concentration,'or-')
minima_plot = plot(min_times, q1(M1/2,TF),'mv', 'MarkerSize',15);
minima_plot.MarkerFaceColor = [1 0.5 0];

title('Midpoint Activator and Repressor Concentrations vs. Time (Synchronous)')
xlabel('Time')
ylabel('Signal Amount')
axis([0 Time/dt 0 1])
ax = gca;
ax.FontSize = 16; 

legend('q1','q2', 'q1 Minima')
saveas(fig1,'P2N1_Initialization_Peaks.pdf')

drawnow

% Calculate phase from time difference between last two minima.
beg_peak = min_times(end-1);
end_peak = min_times(end);
phase_length = round(end_peak-beg_peak);

% Find time corresponding to phase shift
phase_shift_25_time = round(phase_shift_value*phase_length);

% Find concentrations corresponding to 25% phase shift
q1shift25 = q1(M1/2,phase_shift_25_time);
q2shift25 = q2(M1/2,phase_shift_25_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Synchronization Test

% Change to false if you do now want to run the synchronization test
sync_test = true; 

if sync_test == true
    % Reinitialize variables
    Time = 40; %Total Sim Time
    clear q1 q2 q1_del q2_del
    
    % Prepare video
    fig2 = figure(2);
    vidfile = VideoWriter('synchronization_test_p2n1.mp4','MPEG-4');
    open(vidfile);

    % Initialize phases of subpopulations:
    for i=1:(M1/2)
        q1(i,1) = 0; %Initial condition
        q2(i,1) = 0; %Initial condition
    end

    for i=(M1/2)+1:M1
        q1(i,1) = q1shift25; %Initial condition for 25% phase difference.
        q2(i,1) = q2shift25; %Initial condition for 25% phase difference.
    end
    
    % Reinitialize simulation
    for i = 1:M1
        %Setup Tridiagonal Matrix (Diffusion and Decay Terms)
        A(i) = -r_s;
        B(i) = 2.0*r_s+1.0+gamma_s*dt/2.0;
        C(i) = -r_s;
        if (i == 1) % Implement open boundary conditions (set q1,q2 = 0 at boundary)
            A(i) = 0.0;
            B(i) = 1.0;
            C(i) = 0.0;
            F1(i) = 0.0;
            F2(i) = 0.0;
        elseif (i == M1)
            A(i) = 0.0;
            B(i) = 1.0;
            C(i) = 0.0;
            F1(i) = 0.0;
            F2(i) = 0.0;
        end
    end
    
    %Piecewise define production terms so activator production is on for only
    %even stripes and repressor production is on for only odd stripes.
    for i = 1:M1
        if (mod(floor(i/(M1/M)),2) == 0) % CHANGE!
            na0_s(i) = na0p_s;
            na1_s(i) = na1p_s;
            nr0_s(i) = 0.0;
            nr1_s(i) = 0.0;
        else
            na0_s(i) = 0.0;
            na1_s(i) = 0.0;
            nr0_s(i) = nr0p_s;
            nr1_s(i) = nr1p_s;
        end
    end
    
    %Spatial Grid
    for i = 1:M1
        x(i) = (i-1)*dx;
    end
    
    for k = 1:Time/dt
        %Delay Term
        for i = 1:M1
            if (k <= Delay)
                q1_del(i) = 0.0;
                q2_del(i) = 0.0;
            else
                q1_del(i) = q1(i,k-Delay);
                q2_del(i) = q2(i,k-Delay);
            end
        end
        for i = 2:M1-1
        %Compute RHS for both q1 and q2 each time step
        F1(i) = q1(i+1,k)*r_s+q1(i,k)*(1-2.0*r_s-dt*gamma_s/2.0)+q1(i-1,k)*r_s+ dt*(na0_s(i)+na1_s(i)*(q1_del(i))^(n1))/(1+(q1_del(i))^n1+(q2_del(i)/K2_s)^n2) - dt*de_s*(q1(i,k)*(q2_del(i)/Ke_s)^ne)/(1+(q2_del(i)/Ke_s)^ne);
        F2(i) = q2(i+1,k)*r_s+q2(i,k)*(1-2.0*r_s-dt*gamma_s/2.0)+q2(i-1,k)*r_s+ dt*(nr0_s(i)+nr1_s(i)*(q1_del(i))^(n1))/(1+(q1_del(i))^n1) - dt*de_s*(q2(i,k)*(q2_del(i)/Ke_s)^ne)/(1+(q2_del(i)/Ke_s)^ne);
        end
        
        % Solve the tridiagonal matrix to get the next time step values
        q1(:,k+1) = tridiag(B,A,C,F1);
        q2(:,k+1) = tridiag(B,A,C,F2);

        % Plot values and save to video
        if (mod(k,500)==0)
            clf()
            set(fig2,'visible','on');
            hold on
            gca()
            plot(x,q1(:,k+1),'ob-')
            plot(x,q2(:,k+1),'or-')
    
            title('Synchronization Test (Signal Amount vs. Position)')
            xlabel('Position')
            ylabel('Signal Amount')
            axis([0 L 0 0.1*10])
            ax = gca;
            ax.FontSize = 16; 
    
            legend('q1','q2')
            k % PRINTS THE TIME STEP!
            drawnow
    
            F(k) = getframe(gcf);
            writeVideo(vidfile, F(k));
        end        
        
    end

    % Close video file
    close(vidfile)

    % Get q1 concentration for middle left and middle right subpop
    left_concentration = q1(M1/4,:);
    right_concentration = q1(3*M1/4,:);

    % Get normailzed concentrations
    normalized_left_concentration = q1(M1/4,:)/max(q1(M1/4,:));
    normalized_right_concentration = q1(3*M1/4,:)/max(q1(3*M1/4,:));

    % Plot nomalized concentration for middle left and middle right subpop
    % vs. time
    sim_time = 1:Time/dt+1;
    fig3 = figure(3);
    clf()
    set(fig3,'visible','on');
    hold on
    gca()
    plot(sim_time,normalized_left_concentration,'ob-')
    plot(sim_time,normalized_right_concentration,'or-')

    title('P2N1 Signal Production vs. Time')
    xlabel('Time')
    ylabel('Signal Amount')
    axis([0 Time/dt 0 1])
    ax = gca;
    ax.FontSize = 16; 

    legend('a - left half','a - right half')
    saveas(fig3,'P2N1_Concentration_Peaks.pdf')

    % Get minima for left and right side
    left_TF = islocalmin(q1(M1/4,:));
    left_min_times = sim_time(left_TF)
    right_TF = islocalmin(q1(3*M1/4,:));
    right_min_times = sim_time(right_TF)
    true_right_min_times = right_min_times(right_min_times>1000)

    phase_length

    phase_difference_percent = abs(left_min_times(10)-true_right_min_times(10))/phase_length*100;

    fprintf('The percent phase difference for P2N1 is: %.2f\n',phase_difference_percent);

end





