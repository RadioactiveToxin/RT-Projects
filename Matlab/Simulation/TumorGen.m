classdef TumorGen % Tumor Generator
    % this project is based on the research paper article titled:
    % Three-dimensional tumor growth in timevarying chemical fields: a modeling framework and theoretical study
    % https://doi.org/10.1186/s12859-019-2997-9

    properties(GetAccess=public)
        Origin=[5 5 5];
        Dims=[1 1 1];
        N=10;
        alpha=1;
        VO=zeros(10^3,3); % Voxels Origin
        Active=[];
        Mmax=100000;
        M_max=(100000)*ones(10^3,1);
        Boundary=[]; % Boundary Voxels
        ds=1;
        time=4;
        dt=1;
        K=2;
        mu=5;
        lt=zeros(10^3,1);
        nc_t=(0)*ones(10^3,1);
        nn_t=(0)*ones(10^3,1);
    end

    methods
        % Constructor
        function obj=TumorGen(Origin,ds,N,Mmax,alpha,K,mu)
            switch(nargin)
                case 1
                    obj.Origin=Origin;
                case 2
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds; 
                case 3
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds;
                    obj.N=N; obj.lt=zeros(N^3,1); obj.VO=zeros(N^3,3);
                    obj.nc_t=(0)*ones(N^3,1); obj.nn_t=(0)*ones(N^3,1);
                    obj.M_max=(100000)*ones(N^3,1);
                case 4
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds;
                    obj.M_max=Mmax*ones(N^3,1); obj.Mmax=Mmax;
                    obj.N=N; obj.lt=zeros(N^3,1); obj.VO=zeros(N^3,3);
                    obj.nc_t=(0)*ones(N^3,1); obj.nn_t=(0)*ones(N^3,1);
                case 5
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds;
                    obj.M_max=Mmax*ones(N^3,1); obj.Mmax=Mmax; obj.N=N;
                    obj.alpha=alpha; obj.lt=zeros(N^3,1); obj.VO=zeros(N^3,3);
                    obj.nc_t=(0)*ones(N^3,1); obj.nn_t=(0)*ones(N^3,1);
                case 6
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds;
                    obj.M_max=Mmax*ones(N^3,1); obj.Mmax=Mmax; obj.K=K;
                    obj.N=N; obj.alpha=alpha; obj.lt=zeros(N^3,1); obj.VO=zeros(N^3,3);
                    obj.nc_t=(0)*ones(N^3,1); obj.nn_t=(0)*ones(N^3,1);
                case 7
                    obj.Origin=Origin; obj.Dims=[ds ds ds]; obj.ds=ds;
                    obj.M_max=Mmax*ones(N^3,1); obj.Mmax=Mmax; obj.K=K; obj.mu=mu;
                    obj.N=N; obj.alpha=alpha; obj.lt=zeros(N^3,1); obj.VO=zeros(N^3,3);
                    obj.nc_t=(0)*ones(N^3,1); obj.nn_t=(0)*ones(N^3,1);
            end
            for i=1:obj.ds:(obj.N)
                for j=1:obj.ds:(obj.N)
                    for k=1:obj.ds:(obj.N)
                        obj.VO((i+(j-1)*(obj.N)+(k-1)*(obj.N)^2),:)=[i j k];
                        if (i==1)|(i==obj.N)|(j==1)|(j==obj.N)|(k==1)|(k==obj.N)
                            obj.Boundary=[obj.Boundary;i j k];
                        end
                    end
                end
            end
            [~,row]=ismember(obj.Origin,obj.VO,'rows');
            obj.lt(row)=1;
        end

        % Tumor(Voxel) Plot, VoxelPlot(Coordinate,Dimensions,Color,transparency)
        % This is a modefied version of the actual function. 
        % you can find the actual function in the link below:
        % https://www.mathworks.com/matlabcentral/fileexchange/3280-voxel?s_tid=srchtitle_site_search_1_voxel%2520plot
        function []=VoxelPlot(obj,Cordi,D,Color,alpha)
            grid on; view([129.63 21.80]);
            axis([1 (obj.N)+1 1 (obj.N)+1 1 (obj.N)+1]);
            rotate3d on; xlabel('x');
            ylabel('y'); zlabel('z');
            switch(nargin)
                case 1
                    disp('Too few arguements for voxel');
                    return;
                case 2
                    l=1;    % default length of side of voxel is 1
                    Color='b';  % default color of voxel is blue
                case 3
                    Color='b';
                case 4
                    alpha=1;
                case 5
                    % do nothing
                otherwise
                    disp('Too many arguements for voxel');
            end
            x=[Cordi(1)+[0 0 0 0 D(1) D(1) D(1) D(1)];Cordi(2)+[0 0 D(2) D(2) 0 0 D(2) D(2)];Cordi(3)+[0 D(3) 0 D(3) 0 D(3) 0 D(3)]]';
            for n=1:3
                if n==3
                    x=sortrows(x,[n,1]);
                else
                    x=sortrows(x,[n n+1]);
                end
                temp=x(3,:);
                x(3,:)=x(4,:);
                x(4,:)=temp;
                h=patch(x(1:4,1),x(1:4,2),x(1:4,3),Color);
                set(h,'FaceAlpha',alpha);
                temp=x(7,:);
                x(7,:)=x(8,:);
                x(8,:)=temp;
                h=patch(x(5:8,1),x(5:8,2),x(5:8,3),Color);
                set(h,'FaceAlpha',alpha);
            end
        end

        % Algorithm 1: glucose_oxygen Diffusion
        function [T,Tc]=Algorithm1(obj)
            T=zeros((obj.N)^3);
            [m,~]=size(obj.VO);
            P1=9.0009; P2=0.9001; P3=0.0900;
            Pm=sum([P1 P2 P3]);
            for e=1:m
                [~,index]=ismember(obj.VO(e,:),obj.Boundary,'rows');
                if index~=0
                    T(e,e)=1;
                    Tc(e,e)=0;
                else
                    for di=[-obj.ds 0 obj.ds]
                        for dj=[-obj.ds 0 obj.ds]
                            for dk=[-obj.ds 0 obj.ds]
                                [~,index1]=ismember(obj.VO(e,:)+[di dj dk],obj.VO,'rows');
                                if abs(di)+abs(dj)+abs(dk)==0
                                    T(e,e)=90.0090;
                                    Tc(e,e)=0;
                                elseif abs(di)+abs(dj)+abs(dk)==1*obj.ds
                                    T(e,index1)=P1;
                                    Tc(e,e)=(P1/Pm)*100;
                                elseif abs(di)+abs(dj)+abs(dk)==2*obj.ds
                                    T(e,index1)=P2;
                                    Tc(e,e)=(P2/Pm)*100;
                                elseif abs(di)+abs(dj)+abs(dk)==3*obj.ds
                                    T(e,index1)=P3;
                                    Tc(e,e)=(P3/Pm)*100;
                                end
                            end
                        end
                    end
                end
            end
        end
 
        % The random perturbation of matrices T
        function T=Perturbation(obj,T)
            [~,n]=size(T);
            rand_num=randi([1 1+obj.mu],[n 1]);
            T=rand_num.*T;
            for column=1:n
                if sum(T(:,column))~=0
                    T(:,column)=(T(:,column)/sum(T(:,column)))*100;
                end
            end
        end
        
        % Diffusion
        function []=Algorithm2(obj,Tc)
            S0=obj.M_max-obj.nc_t-obj.nn_t;
            Wt=obj.lt+obj.nc_t+obj.nn_t;
            S1=(S0>=zeros((obj.N)^3,1)).*S0;
            S2=((Wt>obj.M_max).*S1)+((Wt<=obj.M_max).*obj.lt);
            obj.lt=S2+Tc*(obj.lt-S2);
        end

        % Voxel State visualizer
        function color=StateColor(obj,Cordi)
            [~,index]=ismember(Cordi,obj.VO,'rows');
            if obj.lt(index)>=(obj.Mmax)
                color=[0 0 1];
            elseif obj.lt(index)>=(obj.Mmax)/2
                color=[0 0.5 1];
            else
                color=[0 1 1];
            end
        end

        % Tumor Visualization
        function []=Generate(obj)
            [T,Tc]=obj.Algorithm1();
            To=obj.Perturbation(T);
            Tgl=obj.Perturbation(T);
            [~,nTc]=size(Tc);
            rand_num=randi([1 1+obj.mu],[nTc 1]);
            Tc=rand_num.*Tc;
            Tc=floor((Tc/sum(unique(Tc)))*100);
            
            % Video Capture
            Video1=VideoWriter("Tumor Simulation Fast.avi","Motion JPEG AVI");
            Video2=VideoWriter("Tumor Simulation Slow.avi","Motion JPEG AVI");
            Video1.FrameRate=60; Video2.FrameRate=60; % Set the frame rate
            Video1.open(); Video2.open();
            for t=0:obj.dt:obj.time
                % The random perturbation of matrices To, Tgl, Tc
                if mod(t,obj.K*obj.dt)==0
                    To=obj.Perturbation(To);
                    Tgl=obj.Perturbation(Tgl);
                    rand_num=randi([1 1+obj.mu],[nTc 1]);
                    Tc=rand_num.*Tc;
                    Tc=floor((Tc/sum(unique(Tc)))*100);
                end
                % glucose_oxygen Diffusion Effect on Voxel quantity
                obj.lt=floor(To*(Tgl*(obj.lt)));
                % Diffusion to other Voxels
                obj.Algorithm2(Tc);
                % Active
                obj.Active=unique([obj.Active;find(obj.lt>1)]);
                for a=1:length(obj.Active)
                    v=obj.VO(obj.Active(a),:);
                    obj.VoxelPlot(v,obj.Dims,obj.StateColor(v),obj.alpha);
                    % Capture the current frame
                    frame1=getframe(gcf);
                    Video1.writeVideo(frame1); % Write the frame to the video
                    pause(1);
                end
                % Capture the current frame
                frame2=getframe(gcf);
                Video2.writeVideo(frame2); % Write the frame to the video
                pause(1);
            end
            Video1.close(); Video2.close();
        end
    end
end
