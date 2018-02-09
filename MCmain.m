%% initialize photons
close all
photon_type = 1;
ang = 20;

% refractive index
    n1 = 1;
    n2 = 1.37;
PhotonNo = [1e4];

mu_a = 0.1; % in cm^-1
mu_s = 100; % in cm^-1
g = 0.9;
dz = 5e-3;
zmax = 1;
nz = floor(zmax/dz)+1;
xmax = 1;
nx = nz;
dx = 2*xmax/nx;
tic
for i = 1:length(PhotonNo)
    for ndif = 1:3
        switch ndif 
            case 1
                ang = 20;
            case 2
                ang = 45;
            case 3
                ang = 70;
        end
        [x,y,z,ux,uy,uz] = initPhoton(photon_type,ang);
        pos = [x,y,z];
        dir = [ux,uy,uz];
        Rsp = ((n1-n2)/(n1+n2))^2;
        Weight = 1 - Rsp;
        
        Absorption = zeros(nx,floor(zmax/dz)+1);
        z_axis = linspace(0,zmax,length(Absorption));
        for n = 1:PhotonNo(i)
            A = MCPhoton(pos,dir,Weight,mu_a,mu_s,g,n1,n2,dz,zmax,dx,xmax);
            Absorption = A + Absorption;
        end
%         if ndif == 1
        figure;
%             semilogy(z_axis, Absorption/mu_a,'-.','linewidth',1.5);ylabel('Fluence [-]');xlabel('Distance (cm)');
%             name = strcat('Photon No.=',num2str(PhotonNo(i)));
        imagesc(Absorption);
        name = strcat('Incident angle: \theta=',num2str(ang));
        title(name);axis off;axis square
        
%         else
%             hold on
%             semilogy(z_axis, Absorption/mu_a,'linewidth',1.5);ylabel('Fluence [-]');xlabel('Distance (cm)');
%             name = strcat('Photon No.=',num2str(PhotonNo(i)));
%             title(name);
%             legend('n_{rel}=1','n_{rel}=1.37');
%         end
        
    end
end
toc