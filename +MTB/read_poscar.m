function obj = read_poscar(obj,filename)
            %read POSCAR
            %Return:
            %a : lattice constants
            %b : Reciprocal vectors
            %
            %b=inv(a')*2*pi % a'*b=2pi
            %
            poscar=readtable(filename);
            a1=table2array(poscar(1,:));
            a2=table2array(poscar(2,:));
            a3=table2array(poscar(3,:));
            a=[a1;a2;a3];
            obj.a=a;
            omega=dot(a1,cross(a2,a3));
            b1=2*pi*cross(a2,a3)/omega;
            b2=2*pi*cross(a3,a1)/omega;
            b3=2*pi*cross(a1,a2)/omega;
            b=[b1;b2;b3];
            obj.b = b;
            obj.atoms = table2array(poscar(4:end,:));
        end