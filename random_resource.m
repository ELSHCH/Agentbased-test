%  The main programm determines the positions of a swarm of individuals
%  on the two dimensional domain of resources, the individuals move randomly
%  such that they are attracted
%  to the center of mass of the whole group and have a mutual force of
%  repulsion for every pair of individuals;
%  the resource domain is updated at every step depending on the positions
%  of grazers 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  L_x, L_y, num_x, num_y are the width and height of the square domain, 
%                numbers of boxes for the domain mesh
%  time_step is the time step (inverse of the frequency of jumping from current position to
%                   the subsequent), for every individual it is the same  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  xbox(num_x),ybox(num_y) the coordinates of the mesh on the domain
%  coord_x(1:n),coord_y(1:n),velocities_x(1:n),velocities_y(1:n) coordinates and velocities of n individuals
%  coord_x_init(1:n),coord_y_init(1:n),velocities_x_init(1:n),velocities_y_init(1:n)
%   the initial coordinates and velocities of n individuals
%  coord_t(1:n,1:2,1:num_steps),velocity_t(1:n,1:2,1:num_steps) are the output matrices 
%              where the trajectories and velocities for all the
%              individuals are saved
%  num_steps   the total number of time steps
%  n           the number of individuals for the 
%  t(1:num_steps) is the time vector  
%  time_span(n) 
%  resource(num_x,num_y),resource_init(num_x,num_y) is the distribution of resources at each time steps 
%  distr_resource(num_x,num_y,num_steps) is the distribution of resources
%                    on the domain for the number of steps = num_steps  
%  position_prey(1:num_x,1:num_y) is the matrix of positions of preys on the domain
%                      it has zeros and ones as entries depending on the absence and presence 
%                      of the prey
%  delta_time is the fraction of time for which the area of available
%      resources increases with rate Delta(Area)/rate_resource
%  Overlap(1:num_x,1:num_y) is the matrix of overlap of preys positions and food matrix at current step 
%  center_m(1:num_steps,1:2), angular_momentum(1:num_steps,1:3),
%                             polariz_vec(1:num_steps,1:2) center of mass, angular momentum and
%                             polarization of group of individuals at every time step
clear;
global L_x L_y num_x num_y time_step n speed_1 speed_2 r_max r_min r_orient sigma tau theta
L_x=100;
L_y=200;
num_x=100;
num_y=200;
r_max=100*L_x/num_x;
r_orient=5*L_x/num_x;
r_min=3*L_x/num_x;
time_step=1; 
xbox=zeros(num_x,1);
ybox=zeros(num_y,1);
coord_x_init=zeros(n,1);
coord_y_init=zeros(n,1);
velocities_x_init=zeros(n,1);
velocities_y_init=zeros(n,1);
num_steps=100;
speed_1=0.5;
speed_2=0.1;
sigma=180/pi; % in degrees
tau=1; % in seconds
theta=360; % 20 in degrees
distr_resource=zeros(num_x,num_y,num_steps);
coord_t=zeros(n,2,num_steps);
velocity_t=zeros(n,2,num_steps);
center_m=zeros(num_steps,2);
angular_momentum=zeros(num_steps,3);
polariz_vec=zeros(num_steps,2);
delta_time=100;
for i=1:num_x
   xbox(i)=i*L_x/num_x;    
end;    
for i=1:num_y
   ybox(i)=i*L_y/num_y;    
end;   
n=50; %the number of preys
c_group=[50 50]; % the initial center of the group
% Initial distribution of the coordinates and velocities of the preys
coord_x_init(1)=c_group(1);
coord_y_init(1)=c_group(2);
i=2;
while i<=n    
    coord_x_init(i)=5*(rand(1)-rand(1))+c_group(1);
    coord_y_init(i)=5*(2+rand(1)-rand(1))+c_group(2);
    if sqrt((coord_x_init(i)-coord_x_init(i-1))^2+...
            (coord_y_init(i)-coord_y_init(i-1))^2)<r_min
    else
        i=i+1;
    end;
end;    
for i=1:n
     velocities_x_init(i)=rand(1)-rand(1);
     velocities_y_init(i)=rand(1)-rand(1);
%     velocities_x_init(i)=1;
%     velocities_y_init(i)=0;
    n_v=sqrt(velocities_x_init(i)^2+velocities_y_init(i)^2); % the norm of the velocities
    if n_v~=0
    velocities_x_init(i)=velocities_x_init(i)/n_v;
    velocities_y_init(i)=velocities_y_init(i)/n_v;
    end;
    coord_t(i,1,1)=coord_x_init(i);
    coord_t(i,2,1)=coord_y_init(i);
    velocity_t(i,1,1)=velocities_x_init(i);
    velocity_t(i,2,1)=velocities_y_init(i);   
end;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial distribution of the resources 
% resources are concentraited around preferential locations
concent=num_x*num_y/20000*10; % number of preferential locations 
rad=5; % radius of resource concentration at each preferential location
resource_init=zeros(num_x,num_y); 
k0=0;
while k0<concent
    ind1=randint(1,1,[1,num_x]);
    ind2=randint(1,1,[1,num_y]);
    resource_init(ind1,ind2)=1;
    k0=k0+1;
    for k1=-rad:1:rad
    for k2=-floor(sqrt(rad^2-k1^2)):1:floor(sqrt(rad^2-k1^2))  
    if  (ind1+k1>=1)&&(ind1+k1<=num_x)&&(ind2+k2<=num_y)&&(ind2+k2>=1)   
      %resource_init(ind1+k1,ind2+k2)=randint(1,1,[0 1]);
      resource_init(ind1+k1,ind2+k2)=1;
    end;    
    end;
    end;
end;
for i1=1:num_x
    for i2=1:num_y  
    distr_resource(i1,i2,1)=resource_init(i1,i2);        
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Start numerical simulation for a given number of time steps
 Overlap1=zeros(num_x,num_y);
 position_prey=zeros(num_x,num_y);
 Prey_pos=zeros(num_x,num_y,num_steps);
 resource=resource_init;
 vel_inside=zeros(num_x,num_y,2);
 vel_in=zeros(num_x,num_y,2);
 t=zeros(num_steps,1);
t(1)=0;
sp=zeros(num_steps,1); % current speed of the group is used as marker for the intersection
                       % with positive resources zone
sp(1)=speed_1;                       
 for i=2:num_steps
     clear position_prey 
     position_prey=zeros(num_x,num_y);
     t(i)=time_step*(i-1);
     [dx,dy]=gradient(resource);
     %position_prey=zeros(num_x,num_y);
     vel_in=vel_inside;
     [coord_x,coord_y,velocities_x,velocities_y,sp(i)]=...
         position_shift(t(i),coord_x_init,coord_y_init,...
         velocities_x_init,velocities_y_init,xbox,ybox,Overlap1,dx,dy);
% Determine the new positions of prey population on the square domain 
    for jj=1:n
    for ii1=1:num_x-1
    for ii2=1:num_y-1
        if (coord_x(jj)<=xbox(ii1+1))&&(coord_x(jj)>xbox(ii1))
        if (coord_y(jj)<=ybox(ii2+1))&&(coord_y(jj)>ybox(ii2))
            position_prey(ii1,ii2)=1;
            vel_inside(ii1,ii2,1)=velocities_x(jj);
            vel_inside(ii1,ii2,2)=velocities_y(jj);
        end;
        end;
    end;    
    end;    
    end;
%     resource_init
%     position_prey
     [resource, Overlap1]=resource_distr_matrix(t(i),delta_time,resource_init,position_prey);
     for i1=1:num_x
    for i2=1:num_y  
    distr_resource(i1,i2,i)=resource(i1,i2); 
    Prey_pos(i1,i2,i)=-1*(position_prey(i1,i2)-1); %places zeros whenever the position is occupied by an individual
    end;
    end;
    [dx,dy]=grad(distr_resource(1:num_x,1:num_y,1));
    for i3=1:n
    coord_t(i3,1,i)=coord_x(i3);
    coord_t(i3,2,i)=coord_y(i3);
    velocity_t(i3,1,i)=velocities_x(i3);
    velocity_t(i3,2,i)=velocities_y(i3);
    end;
    resource_init=resource;
    coord_x_init=coord_x;
    coord_y_init=coord_y;
    velocities_x_init=velocities_x;
    velocities_y_init=velocities_y;
 end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate diffusion rate as number of boxes visited by the entire population per unit time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% unit_time=10*time_step
%
num_probes=floor(num_steps/(10*time_step));
Diffusion=ones(num_x,num_y,num_probes);
for kk=1:num_probes
for is1=1:10
    for ix=1:num_x
    for iy=1:num_y    
      Diffusion(ix,iy,kk)=Diffusion(ix,iy,kk)*Prey_pos(ix,iy,is1+(kk-1)*10);
    end;
    end;
end;
s=nonzeros(Diffusion(:,:,kk));
diffusion_rate(kk)=num_x*num_y-length(s);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate group polarization, center of mass and angular momentum 
[polariz,center_m,angular_momentum]=polarization(velocity_t(1:n,1,1:num_steps),...
          velocity_t(1:n,2,1:num_steps),...
          coord_t(1:n,1,1:num_steps),coord_t(1:n,2,1:num_steps));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  Save the variables in the output .mat file
 save random_resource_output19 distr_resource coord_t t diffusion_rate velocity_t sp ...
     polariz center_m angular_momentum
%  Graphical output
[Xmesh,Ymesh]=meshgrid(xbox,ybox);
mov = avifile('resource_19.avi');
 axis([0 100 0 100]);
     axis(axis);
for i=1:num_steps
   pcolor(Xmesh,Ymesh,distr_resource(1:num_x,1:num_y,i)');
   hold on
   shading flat
%   plot(coord_t(:,1,i),coord_t(:,2,i),'o','MarkerSize',5,'MarkerFaceColor','g');
%    quiver(coord_t(:,1,i),coord_t(:,2,i),...
%          velocity_t(:,1,i),velocity_t(:,2,i));
 %  hold off
   F = getframe(gca);
   mov = addframe(mov,F);
end;
mov = close(mov);