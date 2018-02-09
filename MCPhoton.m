function A = MCPhoton(pos,dir,weight,mu_a,mu_s,g,n1,n2,dz,zmax,dx,xmax)
% This is the main script for photon procedure
% if beamtype == 1    % pencil beam
%     ux = 0;uy = 0;uz = 1;
% end
x = pos(1);y = pos(2);z = pos(3);
ux = dir(1);uy = dir(2);uz = dir(3);
% boundary
z0 = 0; % z1 = infinity
%% free path s
mu_t = mu_a + mu_s;
dead = 0;
iter = 800;
i = 1;
nzmax = floor(zmax/dz)+1;
nxmax = floor(2*xmax/dx);
A = zeros(nxmax,nzmax);

while(dead==0)
    
    % set new s
    s = -log(rand)/(mu_t);
    
    % apply boundary condition
    if uz < 0   % right now only one layer
        db = (z0 - z)/uz;
        if db <= s     % photon hit the boundary
            s = s - db;
            alphai = acos(abs(uz));
            alphat = asin(n2*sin(alphai)/n1);
            if n2 > n1
                crit_ang = asin(n1/n2);
                if alphai > crit_ang
                    R = 1;  % total internal reflectance
                else
                    R = 1/2*(sin(alphai-alphat)^2/sin(alphai+alphat)^2 + tan(alphai-alphat)^2/tan(alphai+alphat)^2);
                end
            else
                R = 1/2*(sin(alphai-alphat)^2/sin(alphai+alphat)^2 + tan(alphai-alphat)^2/tan(alphai+alphat)^2);
            end
            if rand <= R    % internally reflected
                % update position
                x = x + ux*db;
                y = y + uy*db;
                z = z + uz*db;
                uz = -uz;
                
            else
                dead = 1;   % do not record beam going outside tissue
                break;
            end
            
        end
    end
    % update position
    x = x + ux*s;
    y = y + uy*s;
    z = z + uz*s;
    if z < 0
        dead = 1;
        break;
    end
    % Absorption
    dW = (mu_a/mu_t)*weight;
    nz = floor(z/dz)+1;
%     traverse = sqrt(x^2+y^2);
    nx = floor((x+xmax)/dx)+1;
    if nz <= nzmax
        if x <= xmax && x >-xmax
            
            A(nx,nz) = A(nx,nz) + dW;
        end
    end
    weight = weight - dW;
    % photon scattering
    % compute cos_theta
    if g ==0
        cos_theta = 2*rand - 1;
    else
        cos_theta = 1/2/g*(1 + g^2 - ((1-g^2)/(1-g+2*g*rand))^2);
    end
    sin_theta = sqrt(1-cos_theta^2);
    phi = 2*pi*rand;cos_phi = cos(phi); sin_phi = sin(phi);
    if abs(uz)<1-1e-5
        uzz = sqrt(1-uz^2);
        uxnew = sin_theta*(ux*uz*cos_phi - uy*sin_phi)/uzz + ux*cos_theta;
        uynew = sin_theta*(uy*uz*cos_phi + ux*sin_phi)/uzz + uy*cos_theta;
        uznew = -sin_theta*cos_phi*uzz + uz*cos_theta;
    else
        uxnew = sin_theta*cos_phi;
        uynew = sin_theta*sin_phi;
        if uz >0
            uznew = cos_theta;
        else
            uznew = -cos_theta;
        end
    end
    
    
    ux = uxnew;
    uy = uynew;
    uz = uznew;
    % Russian roulette
    if weight <1e-3
        % 1 chance in 10 to survive
        if rand <= 0.1
            weight = 10*weight;
        else
            dead = 1;
            break;
        end
    end

    i = i+1;
end
nz_max = floor(zmax/dz)+1;
nx_max = floor(2*xmax/dz)+1;
A(nx_max:end,nz_max+1:end) = [];
end