classdef Optimize
    properties(GetAccess=public)
        F;
        G;
        C;
        A;
        b;
        X;
        B;
        R;
        Simplex_Matrix;
    end
        
    methods
        % Constructor
        function obj=Optimize()
            prompt={'Objective Function','Constructive Constraints','Sign Constraints'};
            dlgtitle='linear problem'; dims=[2 60;5 60;2 60];
            Input=inputdlg(prompt,dlgtitle,dims);
            f=string(Input{1}); g=string(Input{2}); SignC=split(Input{3},' , ');
            
            % Analysis of f(x)
            if contains(f,"max",'IgnoreCase',true)==1
                f=erase(lower(f),"max");
                f="-("+f+")";
            else
                f=erase(lower(f),"min");
            end

            % Analysis of Sign Constraints
            for i=1:length(SignC)
                variables=erase(SignC{i},' ');
                if contains(variables,'free','IgnoreCase',true)==1
                    variables=erase(variables,'free');
                    sign='free';
                elseif contains(variables,{'<0','0>','0>=','<=0'})==1
                    variables=erase(variables,{'<0','0>','0>=','<=0'});
                    sign='<0';
                else
                    sign='>0';
                end
                for Is=string(split(variables,','))
                    switch sign
                        case 'free'
                            Fs=sym(Is+"%d",[1 2]);
                            Fs="("+string(Fs(2)-Fs(1))+")";
                            f=replace(f,Is,Fs);
                            g=replace(g,Is,Fs);
                        case '<0'
                            f=replace(f,Is,"(-"+replace(Is,"x","y")+")");
                            g=replace(g,Is,"(-"+replace(Is,"x","y")+")");
                    end
                end
            end

            % Analysis of g(x)*
            g=char(g); rhs=[]; rg=[];
            L=size(g); k=0;
            for j=1:L(1)
                gi=g(j,:);
                gi=replace(replace(erase(gi,' '),'<=','<'),'>=','>');
                if contains(gi,{'<','>'})==1
                    if count(gi,{'<','>'})==2
                        k=k+2;
                        syms 's%d' [1 k];
                        syms 'R%d' [k 1];
                        if contains(gi,'<')==1
                            S=strfind(gi,'<');
                            rhs=[rhs;string(gi(S(2)+1:end));string(gi(1:S(1)-1))];
                        else
                            S=strfind(gi,'>');
                            rhs=[rhs;string(gi(1:S(1)-1));string(gi(S(2)+1:end))];  
                        end
                        rg=[rg;string(gi((S(1)+1):(S(2)-1)))+"+"+string(s(k-1))+"+"+string(R(k-1));string(gi((S(1)+1):(S(2)-1)))+"-"+string(s(k))+"+"+string(R(k))];
                    else
                        k=k+1;
                        syms 's%d' [1 k];
                        syms 'R%d' [k 1];
                        if contains(gi,'<')==1
                            S=strfind(gi,'<');
                            if isempty(str2num(string(gi(1:S-1))))==1
                                rhs=[rhs;string(gi(S+1:end))];
                                rg=[rg;string(gi(1:S-1))+"+"+string(s(k))+"+"+string(R(k))];
                            else
                                rhs=[rhs;string(gi(1:S-1))];
                                rg=[rg;string(gi(S+1:end))+"-"+string(s(k))+"+"+string(R(k))];
                            end
                        else
                            S=strfind(gi,'>');
                            if isempty(str2num(string(gi(1:S-1))))==1
                                rhs=[rhs;string(gi(S+1:end))];
                                rg=[rg;string(gi(1:S-1))+"-"+string(s(k))+"+"+string(R(k))];
                            else
                                rhs=[rhs;string(gi(1:S-1))];
                                rg=[rg;string(gi(S+1:end))+"+"+string(s(k))+"+"+string(R(k))];
                            end
                        end
                    end   
                else
                    S=strfind(gi,'=');
                    k=k+1;
                    syms 'R%d' [k 1];
                    if isempty(str2num(string(gi(1:S-1))))==1
                        rg=[rg;string(gi(1:S-1))+"+"+string(R(k))];
                        rhs=[rhs;string(gi(S+1:end))];
                    else
                        rg=[rg;string(gi(S+1:end))+"+"+string(R(k))];
                        rhs=[rhs;string(gi(1:S-1))];
                    end
                end   
            end    
            rg=replace(rg,"s11","s1");
            rg=replace(rg,"R11","R1");
            obj.F=str2sym(f); obj.G=str2sym(rg);
            obj.X=union(symvar(obj.F),symvar(obj.G));
            obj.C=jacobian(obj.F,obj.X); obj.A=jacobian(obj.G,obj.X);
            obj.b=str2double(rhs);
            obj.B=R; obj.R=R;
        end

        function Simplex_Table=ShowTable(obj)
            Simplex_Table=[["NAN" "Z";"Z" 1;obj.B zeros(length(obj.B),1)] [obj.X "RHS";obj.Simplex_Matrix]];
        end

        function M=simplex_pivot(obj,i,s,m)
            obj.Simplex_Matrix(i,:)=obj.Simplex_Matrix(i,:)/obj.Simplex_Matrix(i,s);
            for j=setdiff(1:m,i)
                obj.Simplex_Matrix(j,:)=obj.Simplex_Matrix(j,:)-(obj.Simplex_Matrix(j,s)*obj.Simplex_Matrix(i,:));
            end
            M=obj.Simplex_Matrix;
        end

        function []=Simplex(obj)
            [m,n]=size(obj.A);
            obj.Simplex_Matrix=[sum([-ones(1,m) zeros(1,n-m+1);obj.A obj.b]);obj.A obj.b];
            [m,n]=size(obj.Simplex_Matrix);
            ShowTable(obj)

            % first Phaze
            flag="not optimized";
            while flag=="not optimized"
                % Bland -> min
                s=min(find(obj.Simplex_Matrix(1,1:end-1)==max(obj.Simplex_Matrix(1,1:end-1))));
                if obj.Simplex_Matrix(1,s)<=0
                    disp("Lp1 is optimized");
                    flag="optimized";
                elseif 0<obj.Simplex_Matrix(1,s) && all(obj.Simplex_Matrix(2:end,s)<=0)
                    disp("Lp1 is infinite");
                    flag="infinite";
                else
                    % Bland -> min
                    MIN=obj.Simplex_Matrix(2:end,end)./obj.Simplex_Matrix(2:end,s);
                    i=min(find(MIN==min(MIN(MIN>=0))))+1;
                    obj.Simplex_Matrix=simplex_pivot(obj,i,s,m);
                    obj.B(i-1)=obj.X(s);
                    ShowTable(obj)
                end
            end

            if flag=="optimized"
                % remove Ri
                if obj.Simplex_Matrix(1,end)>0
                    disp("Lp is not feasible");
                else
                    if length(intersect(obj.R,obj.B))~=0
                        L=intersect(obj.R,obj.B);
                        for e=L
                            q=find(obj.B==intersect(obj.B,e));
                            if all(obj.Simplex_Matrix(q+1,m:end-1)==0)
                                obj.Simplex_Matrix(q+1,:)=[];
                                obj.B(q)=[];
                                [m,n]=size(obj.Simplex_Matrix);
                                ShowTable(obj)
                            else
                                p=m-1+find(obj.Simplex_Matrix(q+1,m:end-1)~=0,1,'first');
                                obj.Simplex_Matrix=simplex_pivot(obj,q+1,p,m);
                                obj.B(q)=obj.X(p);
                                ShowTable(obj)
                            end
                        end
                    end
                    obj.Simplex_Matrix(:,1:length(obj.R))=[];
                    obj.C(1:length(obj.R))=[];
                    obj.X(1:length(obj.R))=[];
                    [m,n]=size(obj.Simplex_Matrix);

                    % C
                    obj.Simplex_Matrix(1,:)=[-(obj.C) 0];
                    ShowTable(obj)
                    for k=1:length(obj.B)
                        u=find(obj.X==intersect(obj.X,obj.B(k)));
                        v=find(obj.B==intersect(obj.B,obj.B(k)));
                        obj.Simplex_Matrix(1,:)=obj.Simplex_Matrix(1,:)-obj.Simplex_Matrix(1,u)*obj.Simplex_Matrix(v+1,:);
                    end
                    ShowTable(obj)

                    % second Phaze
                    flag="not optimized";
                    while flag=="not optimized"
                        % Bland -> min
                        s=min(find(obj.Simplex_Matrix(1,1:end-1)==max(obj.Simplex_Matrix(1,1:end-1))));
                        if obj.Simplex_Matrix(1,s)<=0
                            disp("Lp is optimized");
                            flag="optimized";
                        elseif 0<obj.Simplex_Matrix(1,s) && all(obj.Simplex_Matrix(2:end,s)<=0)
                            disp("Lp is infinite");
                            flag="infinite";
                        else
                            % Bland -> min
                            MIN=obj.Simplex_Matrix(2:end,end)./obj.Simplex_Matrix(2:end,s);
                            i=min(find(MIN==min(MIN(MIN>=0))))+1;
                            obj.Simplex_Matrix=simplex_pivot(obj,i,s,m);
                            obj.B(i-1)=obj.X(s);
                            ShowTable(obj)
                        end
                    end
                end     
            end
        end
    end
end

% Examples
% min x1 + 3*x2
% 2*x1 +x2 >= 2
% x1 + 3*x2 >= 3
% x1,x2>0

% min -x1 + 2*x3
% 2*x1 - 4*x2 +x3 <= 4
% x1-3*x2 >=2
% x1,x2,x3>0

% min -4*x1 +2*x2 -6*x3 + 3*x4
% x1 - 2*x2 +3*x3 +x4 = 3
% 2*x1 +2*x2 -x3 -2*x4 = 9
% 3*x1 +2*x3 -x4 = 6
% x1,x2,x3,x4>0

