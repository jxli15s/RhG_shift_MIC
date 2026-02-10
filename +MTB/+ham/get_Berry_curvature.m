function [Omega_k,KX,KY] = get_Berry_curvature(bandindex,Unk,KX,KY,plottap)
%Berry_curvature: This function cal the berrycurvature with the loop method, only work for isolated bands     
               % knum: size of the kmesh
               % bandindex: 1 or 2 or 3... for one bands only 
               % Unk: Norbit*Nband*knum*knum , the eigenstate
            knum=size(Unk,3);
            Omega_k = zeros(knum - 1, knum - 1);            
            Unk=Unk(:,bandindex,:,:);          
            for i = 1:knum - 1
                for j = 1:knum - 1
                    Omega_k(i, j) = small_loop(i, j, Unk);
                end
            end    
            Omega_k=Omega_k./2./pi;
            function Omega = small_loop(i, j, Unk)
                Omega1 = Link(Unk(:,:,i, j), Unk(:,:,i + 1, j)) * Link(Unk(:,:,i + 1, j), Unk(:,:,i + 1, j + 1));
                Omega2 = Link(Unk(:,:,i, j + 1), Unk(:,:,i + 1, j + 1)) * Link(Unk(:,:,i, j), Unk(:,:,i, j + 1));
                Omega = angle(Omega1 / Omega2);
                function U = Link(u1, u2)
                    U = det(u1' * u2) / norm(det(u1' * u2));
                end

            end
           % Plot the Berry Curvature
            KX=KX(1:knum-1,1:knum-1);KY=KY(1:knum-1,1:knum-1);
            % if plottap==1
            % figure('Color','white')
            % contourf(KX,KY,Omega_k)
            % axis equal
            % colorbar
            % % colormap(flipud(gray))
            % title('$\Omega_{k}$','Interpreter','latex')
            % set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
            % end

            if plottap==1
            figure('Color','white')
            % contourf(KX,KY,Omega_k)
            s=pcolor(KX,KY,Omega_k);
            s.EdgeColor='none';
            axis equal
            colorbar
            % colormap(flipud(gray))
            title('$\Omega_{k}$','Interpreter','latex')
            set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
            end
            
        end