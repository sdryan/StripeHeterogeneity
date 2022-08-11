clear
%tridiag(a,b,c,f)

% Values for Crank-Nicolson Grid Initialization
D = 4080;  %Diffusion Strength
tau = 7.5;
L = 4000/(sqrt(tau*D)); %Length of domain
M1 = 7140; %Number of spatial grid points
dx = L/M1; %Spatial Step
dt = .0001;  %Time Step
Time = 20; %Total Sim Time
Track = 0; %Track images
Delay = 10000; %Delay Time (Must be Integer indicating time step delay)
r = D*dt/(2.0*dx*dx); %CN Parameter

%Parameters
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

% Initialize heterogeneous stripe widths
left_stripes = 84;
right_stripes = 14;
left_stripe_types = [];
right_stripe_types = [];
spatial_points_half = M1/2;

% Initialize left stripes
grad_left = spatial_points_half/(left_stripes*(left_stripes+1)/2);
for i = 1:left_stripes
    points_stripe = i*grad_left;
    if mod(i,2) == 1
        stripe_type = 1;
    else
        stripe_type = 2;
    end
    for k = 1:points_stripe
        left_stripe_types(end+1) = stripe_type;
    end
end

% Initialize right stripes
grad_right = spatial_points_half/(right_stripes*(right_stripes+1)/2);
for i = 1:right_stripes
    points_stripe = i*grad_right;
    if mod(i,2) == 1
        stripe_type = 1;
    else
        stripe_type = 2;
    end
    for k = 1:points_stripe
        right_stripe_types(end+1) = stripe_type;
    end
end

% Assemble consortium
stripe_types = [left_stripe_types right_stripe_types];

%Piecewise define production terms so activator production is on for only
%even stripes and repressor production is on for only odd stripes.
for i = 1:M1
    if (stripe_types(i) == 2) % CHANGE!
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

% Prepare video
fig2 = figure(2);
vidfile = VideoWriter('heterogeneity_test_p2n1.mp4','MPEG-4');
open(vidfile)

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
