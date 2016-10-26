function TFPatterns_2D_alpha
close all; clear all

% 
N     = 10;                      % Number of radial mesh points
M     = 40;                      % Number of angular mesh points

t_end   = 50;                      % Simulation time
r_end   = 1;                        % Radial length of domain
del_r   = r_end/(N);            % Radial mesh spacing
del_a   = 2*pi/(M);            % Angular mesh spacing
del_t   = 0.1;                     % Time mesh spacing
    
r       = (del_r/2:del_r:r_end)';            % Radial mesh, r=0 is special
a       = (0:del_a:(2*pi-del_a))';     % Angular mesh (theta = 2pi = 0)
t       = 0:del_t:t_end;            % Time mesh

%diffusion coefficient
D1 = 0.02;
D2 = 0.001;

%reaciton parameters
alpha = 1;
beta = 1;
mu = 0.1;
delta = 0.01;

%gene location for g1
theta = pi/2;
radius = 0.5;
[x,y] = pol2cart(theta,radius);
sprintf('g_1 @ (%g,%g)',x,y)
g1  =   N*find(min(abs(a-radius)))+find(min(abs(a-radius)));
%gene location for g2
theta = pi/8;
radius = 0.9;
[x,y] = pol2cart(theta,radius);
sprintf('g_2 @ (%g,%g)',x,y)
g2  =   N*find(min(abs(a-radius)))+find(min(abs(a-radius)));

%g1 =

p1=ones(N*M+1,1);
p2=ones(N*M+1,1);
p3=ones(N*M+1,1);

p1 = 10*rand(N*M+1,1);
p2 = 10*rand(N*M+1,1);
p3  = 10*rand(N*M+1,1);
%p2(40:60) = 4;
%p1(1:end/2) = 5;
%p2(end/4:3*end/4) = 1.1;

%==========================================================================
% Crank Nicholson Matrices
%==========================================================================
CritCoords = [N*(M/4-1)+1 N*(2*M/4-1)+1 N*(3*M/4-1)+1 N*(M-1)+1];
% Bits and pieces for 2D finite diff polar scheme
u1      = 0.5.*del_t/(del_r.*del_r);
u2      = 0.5.*del_t./(2.*r.*del_r);
u3      = 0.5*del_t./(r.*r.*del_a.*del_a);

F       = zeros(N*M); % 'Forward' matrix
B       = zeros(N*M); % 'Backward' matrix

%NODES
I_C       = zeros(N*M,1);       % Self
I_N       = zeros(N*M-1,1);     % North
I_S       = zeros(N*M-1,1);     % South NEED -1 for each index

for i=1:M
    I_C(N*(i-1)+(1:N))        = (-2*u1-2.*u3);        % All follow same program
    I_N(N*(i-1)+(1:(N-1)))    = (u1-u2(1:(N-1)));     % Excludes top node (no flux)
    I_N(N*(i-1)+1)            = 2*u1;
    I_S(N*(i-1)+(2:(N-1))-1)  = (u1+u2(2:(N-1)));     % Excludes both top an bottom
    I_S(N*i-1)                = 2*u1;                 % No flux at top
    if i~=M; I_S(N*i+1-1) = 0;  end               % Ignore bottom (for now)  
end

I_N(CritCoords) = u1;
I_EW_Main       = zeros(N*(M-1),1); % Nodes East and West, Main
I_EW_Sub       = ones(N,1).*u3;     % Nodes East and West, Sub
for i=1:M-1
    I_EW_Main(N*(i-1)+(1:N)) = u3;
end

%==========================================================================
% for p1
% Make the beasts
F1(1:N*M,1:N*M) = diag(1+D1.*I_C,0)+D1.*( diag(I_S,-1)+diag(I_N,+1)+diag(I_EW_Main,N)+diag(I_EW_Main,-N)+diag(I_EW_Sub,+N*(M-1))+diag(I_EW_Sub,-N*(M-1)));
B1(1:N*M,1:N*M) = diag(1-D1.*I_C,0)+D1.*(-diag(I_S,-1)-diag(I_N,+1)-diag(I_EW_Main,N)-diag(I_EW_Main,-N)-diag(I_EW_Sub,+N*(M-1))-diag(I_EW_Sub,-N*(M-1)));
%clear A B C D E

% Dealing with the center r = 0 'bottom'
F1(N*M+1,CritCoords)       = u1.*D1;        % center axis getting stuff from bottom (M nodes)
B1(N*M+1,CritCoords)       = -u1.*D1;
F1(CritCoords,N*M+1)       = u1.*D1; %bottom getting stuff from center
B1(CritCoords,N*M+1)       = -u1.*D1;
F1(N*M+1,N*M+1)           = 1-4.*u1.*D1;    % center axis losing stuff (M nodes)
B1(N*M+1,N*M+1)            = 1+4*u1.*D1;
B1 = sparse(B1); F1 = sparse(F1);
%==========================================================================
F2(1:N*M,1:N*M) = diag(1+D2.*I_C,0)+D2.*( diag(I_S,-1)+diag(I_N,+1)+diag(I_EW_Main,N)+diag(I_EW_Main,-N)+diag(I_EW_Sub,+N*(M-1))+diag(I_EW_Sub,-N*(M-1)));
B2(1:N*M,1:N*M) = diag(1-D2.*I_C,0)+D2.*(-diag(I_S,-1)-diag(I_N,+1)-diag(I_EW_Main,N)-diag(I_EW_Main,-N)-diag(I_EW_Sub,+N*(M-1))-diag(I_EW_Sub,-N*(M-1)));
%clear A B C D E

% Dealing with the center r = 0 'bottom'
F2(N*M+1,CritCoords)       = u1.*D2;        % center axis getting stuff from bottom (M nodes)
B2(N*M+1,CritCoords)       = -u1.*D2;
F2(CritCoords,N*M+1)       = u1.*D2; %bottom getting stuff from center
B2(CritCoords,N*M+1)       = -u1.*D2;
F2(N*M+1,N*M+1)           = 1-4.*u1.*D2;    % center axis losing stuff (M nodes)
B2(N*M+1,N*M+1)            = 1+4*u1.*D2;
B2 = sparse(B2); F2 = sparse(F2);
%==========================================================================
 



p1_plot = zeros(N+1,M);
p2_plot = zeros(N+1,M);
p3_plot = zeros(N+1,M);


r = [eps; r];
[a,r] = meshgrid(a,r);
[x,y] = pol2cart(a,r);

for n=1:length(t)
    
    %p1(end) = 3;
    p1 =B1\(F1*p1);
    p2 =B2\(F2*p2);
    p3 =B1\(F1*p3);
    %React
    p1  = p1+del_t.*(-delta.*p1+beta.*p2.*p3-alpha.*p1);
    p2  = p2+del_t.*(-delta*p2);
    p3  = p3+del_t.*(-delta.*p3-beta.*p2.*p3+alpha.*p1);
    %p1  =  p1+del_t.*(125.*(0.1-p1+p1.*p1.*p2));
    %p2  =  p2+del_t.*(125.*(1-p1.*p1.*p2));
   
    for i=1:M
        p1(N*i) = p1(N*i)+del_t.*mu.*p3(g1);
        p2(N*i) = p2(N*i)+del_t.*mu.*p3(g2);
    end
    
    if mod(20,10) ==0
    for i=1:M
        p1_plot(2:N+1,i) = p1(N*(i-1)+(1:N));
         p2_plot(2:N+1,i) = p1(N*(i-1)+(1:N));
        p3_plot(2:N+1,i) = p3(N*(i-1)+(1:N));
    end
    p2_plot(1,:) = p1(end); 
    p3_plot(1,:) = p3(end); 
    %p1_plot(1,:) = p1(end); 
    figure(1); %subplot(3,2,1); surf(x, y, p1_plot); zlim([0 3]); subplot(3,2,2); surf(x, y, p1_plot);view(2);  pause(0); drawnow; 
    subplot(2,2,1); surf(x, y, p3_plot); title('Activated TF p1'); zlim([0 round(1.2*max(p3))]); subplot(2,2,2); surf(x, y, p3_plot);view(2);pause(0); drawnow;
    subplot(2,2,3); surf(x, y, p2_plot); title('Inactivator P2'); zlim([0 round(1.2*max(p2))]); subplot(2,2,4); surf(x, y, p2_plot);view(2);pause(0); drawnow;
    % figure(2); surf(x, y, p2_plot);view(2);  pause(0); drawnow; 
    end
end




end