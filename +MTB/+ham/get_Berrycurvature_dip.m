function [Omega_dk,Kx,Ky,Kz]=get_Berrycurvature_dip(obj,Kx,Ky,Kz,Enk,Unk,bandindex,delta,plottap)
%{2023-11-11. cal the berry curvature with quantum geo tensor
     %Kx,Ky: kmesh; obj: the continuum model objetc; 
     %Enk: (eigenvalue) knum*knum*Norbits; 
     %Unk: (eigenvector) Norbits*Nband*knum*knum
     %bandindex: eg: 1:3; delta: small term to avoid the gapless points
     %delta: small term to avoid gapless point
%}                   
            dimH=size(Unk,1);
            % bandindex=dimH+1-bandindex(end):dimH+1-bandindex(1);
            knum=size(Unk,3);
            nband=size(bandindex,2); 
            dkx=obj.b(1,:)/knum;
            dky=obj.b(2,:)/knum;
            dS=norm(dkx)*norm(dky);
            % dS=norm(cross(dkx,dky));
            qm=@MTB.ham.get_quan_metric;            
            dimG=[knum,knum,nband];
            itotal=prod(dimG);
            omega=zeros(1,itotal);
            band_start=bandindex(1);
            % omega_d1=zeros(1,itotal);
            % omega_d2=zeros(1,itotal);
            parfor l=1:itotal
                [i,j,band]=ind2sub(dimG,l);
                fprintf("Processing on %d %d %d of %d KPOINTS\n", i, j, band, itotal);
                band=band_start+band-1;
                k=[Kx(i,j),Ky(i,j),Kz(i,j)]+dkx;
                omega_d1=-2.*imag(qm(obj,k,dkx,dky,Unk(:,:,i,j),reshape(Enk(i,j,:),[dimH,1]),band,delta));
                k=[Kx(i,j),Ky(i,j),Kz(i,j)];
                omega_d2=-2.*imag(qm(obj,k,dkx,dky,Unk(:,:,i,j),reshape(Enk(i,j,:),[dimH,1]),band,delta));
                omega(l)=(omega_d1-omega_d2)./norm(dkx);
            end
            omega=reshape(omega,dimG);
            Omega_dk=omega.*dS;%./2./pi;
            % Omega_k=zeros(knum,knum,size(bandindex,2));
            % for i=1:size(bandindex,2)
            %     Omega_k(:,:,i)=Omega(:,:,n-i);
            % end
            if plottap==1

            figure('Color','white')
            for i=1:nband
                subplot(1,nband,i)
                contourf(Kx,Ky,Omega_dk(:,:,i))
                axis equal
                colorbar
                % colormap(flipud(gray))
                title('$\Omega_{dk}$','Interpreter','latex')
                set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
            end
            end
        end
