function [Omega_k,Kx,Ky]=Berrycurvature_cop(obj,Kx,Ky,Enk,Unk,bandindex,delta,plottap)
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
            [G1,G2]=obj.reciprocal_vectors();
            dkx=G1./knum;
            dky=G2./knum;
            dS=norm(dkx)*norm(dky);
            qm=@continuum.Topology.quan_metric;            
            dimG=[knum,knum,nband];
            itotal=prod(dimG);
            omega=zeros(1,itotal);
            for l=1:itotal
                        [i,j,band]=ind2sub(dimG,l);
                        % band=dimH+1-band;
                        k=[Kx(i,j),Ky(i,j)];                                                                                 
                        omega(l)=-2.*imag(qm(obj,k,dkx,dky,Unk(:,:,i,j),reshape(Enk(i,j,:),[dimH,1]),band,delta));                    
            end      
            omega=reshape(omega,dimG);
            Omega_k=omega.*dS./2./pi;
            % Omega_k=zeros(knum,knum,size(bandindex,2));
            % for i=1:size(bandindex,2)
            %     Omega_k(:,:,i)=Omega(:,:,n-i);
            % end
            if plottap==1

            figure('Color','white')
            for i=1:nband
                subplot(1,nband,i)
                contourf(Kx,Ky,Omega_k(:,:,i))
                axis equal
                colorbar
                % colormap(flipud(gray))
                title('$\Omega_{k}$','Interpreter','latex')
                set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
            end
            end
        end
