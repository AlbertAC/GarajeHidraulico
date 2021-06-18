function ProyectoFinal
close all;
loaddata;
inithome;
n=30;      
    
    K_p = uipanel(fig_1,... %creacion del panel
    'units','pixels',...
    'Position',[20 150 500-120 120],...
    'Title','Cinemàtica Directa','FontSize',11,'BackgroundColor',[.2 .5 .7]);


    LD = 80; % Left, used to set the GUI.

    HT = 18; % Height

    BT = 70; % Bottom
    
    q1_slider = uicontrol(K_p,'style','slider',...
        'Max',6,'Min',0,'Value',0,...
        'SliderStep',[0.01 0.5],...
        'callback',@q1_slider_button_press,...
        'Position',[LD BT 120 HT]);

    q1_min = uicontrol(K_p,'style','text',...
        'String','0',...
        'Position',[LD-30 BT+1 25 HT-4]);

    q1_max = uicontrol(K_p,'style','text',...
        'String','10',...
        'Position',[LD+125 BT+1 25 HT-4]);
    %% loaddata
    function loaddata  
    % La base
        [VCB, FCB, m, stltitle] = stlread('cubo.stl');
         VCB = [VCB*2 ones(size(VCB,1),1)]; %convierte a coordenadas homogeneas
         VCB = VCB*transl(0,0,0);
    % Cabina
        [VCC, FCC, m, stltitle] = stlread('cab.stl');
         VCC = [VCC*2 ones(size(VCC,1),1)]; %convierte a coordenadas homogeneas
         VCC = VCC*transl(0,0,-6);
    % Eslabon 1-1
        [VE11, FE11, m, stltitle] = stlread('base.stl');
         VE11 = [VE11*2 ones(size(VE11,1),1)]; %convierte a coordenadas homogeneas
         VE11 = VE11*transl(7,-7,0);
    % Eslabon 1-2 
        [VE12, FE12, m, stltitle] = stlread('base.stl');
         VE12 = [VE12*2 ones(size(VE12,1),1)]; %convierte a coordenadas homogeneas
         VE12 = VE12*transl(7,-7,0);
    % Eslabon 1-3 
        [VE13, FE13, m, stltitle] = stlread('base.stl');
         VE13 = [VE13*2 ones(size(VE13,1),1)]; %convierte a coordenadas homogeneas
         VE13 = VE13*transl(7,-7,0);
    % Eslabon 1-4    
        [VE14, FE14, m, stltitle] = stlread('base.stl');
         VE14 = [VE14*2 ones(size(VE14,1),1)]; %convierte a coordenadas homogeneas
         VE14 = VE14*transl(7,-7,0);
    % Eslabon 2-1
        [VE21, FE21, m, stltitle] = stlread('eslabon.stl');
         VE21 = [VE21*2 ones(size(VE21,1),1)]; %convierte a coordenadas homogeneas
         VE21 = VE21*transl(7,-7,-6.5);
    % Eslabon 2-2 
        [VE22, FE22, m, stltitle] = stlread('eslabon.stl');
         VE22 = [VE22*2 ones(size(VE22,1),1)]; %convierte a coordenadas homogeneas
         VE22 = VE22*transl(7,-7,-6.5);
    % Eslabon 2-3 
        [VE23, FE23, m, stltitle] = stlread('eslabon.stl');
         VE23 = [VE23*2 ones(size(VE23,1),1)]; %convierte a coordenadas homogeneas
         VE23 = VE23*transl(7,-7,-6.5);
    % Eslabon 2-4    
        [VE24, FE24, m, stltitle] = stlread('eslabon.stl');
         VE24 = [VE24*2 ones(size(VE24,1),1)]; %convierte a coordenadas homogeneas
         VE24 = VE24*transl(7,-7,-6.5);
    % Car   
        [VECR, FECR, m, stltitle] = stlread('auto.stl');
         VECR = [VECR*(1.6) ones(size(VECR,1),1)]; %convierte a coordenadas homogeneas
         VECR = VECR*transl(0,0,0)*Trotz(90)*transl(0,-0.3,0.5);
         
    % House   
        [VHS, FHS, m, stltitle] = stlread('base.stl');
         VHS = [VHS*(0.5) ones(size(VHS,1),1)]; %convierte a coordenadas homogeneas
         VHS = VHS*transl(0,0,0);
         
         L1 = 2.4449;
         setappdata(0,'L1',L1);
         
         L4 = 1.844;
         setappdata(0,'L4',L4);
         
         offset = -6;
         setappdata(0,'offset',offset);
         
         setappdata(0,'VCB',VCB);
         setappdata(0,'FCB',FCB);
         setappdata(0,'VCC',VCC);
         setappdata(0,'FCC',FCC);
         setappdata(0,'VE11',VE11);
         setappdata(0,'FE11',FE11);
         setappdata(0,'VE12',VE12);
         setappdata(0,'FE12',FE12);
         setappdata(0,'VE13',VE13);
         setappdata(0,'FE13',FE13);
         setappdata(0,'VE14',VE14);
         setappdata(0,'FE14',FE14);
         setappdata(0,'VE21',VE21);
         setappdata(0,'FE21',FE21);
         setappdata(0,'VE22',VE22);
         setappdata(0,'FE22',FE22);
         setappdata(0,'VE23',VE23);
         setappdata(0,'FE23',FE23);
         setappdata(0,'VE24',VE24);
         setappdata(0,'FE24',FE24);
         setappdata(0,'VECR',VECR);
         setappdata(0,'FECR',FECR);
         setappdata(0,'VHS',VHS);
         setappdata(0,'FHS',FHS);
         
    end%loaddata

    %% initHome
    function inithome
        VCB = getappdata(0, 'VCB');
        FCB = getappdata(0, 'FCB');
        VCC = getappdata(0, 'VCC');
        FCC = getappdata(0, 'FCC');
        VE11 = getappdata(0, 'VE11');
        FE11 = getappdata(0, 'FE11');
        VE12 = getappdata(0, 'VE12');
        FE12 = getappdata(0, 'FE12');
        VE13 = getappdata(0, 'VE13');
        FE13 = getappdata(0, 'FE13');
        VE14 = getappdata(0, 'VE14');
        FE14 = getappdata(0, 'FE14');
        VE21 = getappdata(0, 'VE21');
        FE21 = getappdata(0, 'FE21');
        VE22 = getappdata(0, 'VE22');
        FE22 = getappdata(0, 'FE22');
        VE23 = getappdata(0, 'VE23');
        FE23 = getappdata(0, 'FE23');
        VE24 = getappdata(0, 'VE24');
        FE24 = getappdata(0, 'FE24');
        VECR = getappdata(0, 'VECR');
        FECR = getappdata(0, 'FECR');
        VHS = getappdata(0, 'VHS');
        FHS = getappdata(0, 'FHS');
        
        L1= getappdata (0,'L1'); 
        L4= getappdata (0,'L4');
        offset = getappdata(0, 'offset');
        
        
        qOld = [0 0 0]; 
        setappdata(0,'qOld',qOld);
        
        set(0,'units','pixels');
        dim = get(0,'screensize');
%% FIGURAS
%         fig_b = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','BASE 2',...
%             'NumberTitle','off');
%         hold on
%         title('BASE');
%         ejes
%         patch('faces', FCB, 'vertices' ,VCB(:,1:3),'EdgeColor','none','facec', 1/255*[255, 160, 122]);
%         
%        fig_cab = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','CABINA',...
%             'NumberTitle','off');
%        hold on
%        title('BASE');
%        ejes
%        patch('faces', FCC, 'vertices' ,VCC(:,1:3),'EdgeColor','none','facec', 1/255*[255, 160, 122]); 
%          
%         fig_e1 = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','1ER ESLABON',...
%             'NumberTitle','off');
%         hold on
%         title('1ER ESLABON');
%         ejes
%         patch('faces', FE11, 'vertices' ,VE11(:,1:3),'EdgeColor','none','facec', 1/255*[119, 136, 153]); 
%          
%           fig_e2 = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','2do ESLABON',...
%             'NumberTitle','off');
%           hold on
%           title('2do ESLABON');
%           ejes
%           patch('faces', FE21, 'vertices' ,VE21(:,1:3),'EdgeColor','none', 'facec', 1/255*[173, 255, 47]);
%     
%         fig_car = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','1ER ESLABON',...
%             'NumberTitle','off');
%         hold on
%         title('CAR');
%         ejes
%         patch('faces', FECR, 'vertices' ,VECR(:,1:3),'EdgeColor','none','facec', 1/255*[119, 136, 153]); 
%     
%         fig_house = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
%             'Name','1ER ESLABON',...
%             'NumberTitle','off');
%         hold on
%         title('CAR');
%         ejes
%         patch('faces', FHS, 'vertices' ,VHS(:,1:3),'EdgeColor','none','facec', 1/255*[119, 136, 153]); 
    
% Imagen ROBOT completo
       fig_1 = figure('doublebuffer','on','Position',[0,35,dim(3)-200,dim(4)-110],...
            'Name','EXAMEN 2DA FASE',...
            'NumberTitle','off');
        hold on
        title('ROBOT');
        ejes% llena el grafico
        view(158,65); 
         
%% CUBO
         LCB = patch('faces', FCB, 'vertices' ,VCB(:,1:3)); 
         set(LCB, 'facec', 1/255*[255, 160, 122]);
         set(LCB, 'EdgeColor','none');
        
%% CABINA
         Vwork = VCC;
         A_CC = transl(0,0,0);
         Vwork = Vwork * A_CC;
         LCC = patch('faces', FCC, 'vertices' ,Vwork(:,1:3)); 
         set(LCC, 'facec', 1/255*[255, 160, 122]);
         set(LCC, 'EdgeColor','none');
%% ESLABONES 1
        % ESLABON 11
         Vwork = VE11;
         A_11 = transl(-7,-7,0);
         Vwork = Vwork * A_11;
         LE11 = patch('faces', FE11, 'vertices' ,Vwork(:,1:3)); 
         set(LE11, 'facec', 1/255*[119, 136, 153]);
         set(LE11, 'EdgeColor','none'); 
        % ESLABON 12
         Vwork = VE12;
         A_12 = transl(7,-7,0);
         Vwork = Vwork * A_12;
         LE12 = patch('faces', FE12, 'vertices' ,Vwork(:,1:3)); 
         set(LE12, 'facec', 1/255*[119, 136, 153]);
         set(LE12, 'EdgeColor','none'); 
        % ESLABON 13
         Vwork = VE13;
         A_13 = transl(-7,7,0);
         Vwork = Vwork * A_13;
         LE13 = patch('faces', FE13, 'vertices' ,Vwork(:,1:3)); 
         set(LE13, 'facec', 1/255*[119, 136, 153]);
         set(LE13, 'EdgeColor','none'); 
        % ESLABON 14
         Vwork = VE14;
         A_14 = transl(7,7,0);
         Vwork = Vwork * A_14;
         LE14 = patch('faces', FE14, 'vertices' ,Vwork(:,1:3)); 
         set(LE14, 'facec', 1/255*[119, 136, 153]);
         set(LE14, 'EdgeColor','none'); 
%% ESLABONES 2
        % ESLABON 21
         Vwork = VE21;
         A_21 = transl(0,0,0.1);
         A_M1 = A_11 * A_21;
         Vwork = Vwork * A_M1;
         LE21 = patch('faces', FE21, 'vertices' ,Vwork(:,1:3)); 
         set(LE21, 'facec', 1/255*[119, 0, 153]);
         set(LE21, 'EdgeColor','none'); 
         [LXe21,LYe21,LZe21]=ejesF(A_M1,2);
        % ESLABON 22
         Vwork = VE22;
         A_22 = transl(0,0,0.1);
         A_M2 = A_12 * A_22;
         Vwork = Vwork * A_M2;
         LE22 = patch('faces', FE22, 'vertices' ,Vwork(:,1:3)); 
         set(LE22, 'facec', 1/255*[119, 0, 153]);
         set(LE22, 'EdgeColor','none');
         [LXe22,LYe22,LZe22]=ejesF(A_M2,2);
        % ESLABON 23
         Vwork = VE23;
         A_23 = transl(0,0,0.1);
         A_M3 = A_13 * A_23;
         Vwork = Vwork * A_M3;
         LE23 = patch('faces', FE23, 'vertices' ,Vwork(:,1:3)); 
         set(LE23, 'facec', 1/255*[119, 0, 153]);
         set(LE23, 'EdgeColor','none'); 
         [LXe23,LYe23,LZe23]=ejesF(A_M3,2);
        % ESLABON 24
         Vwork = VE24;
         A_24 = transl(0,0,0.1);
         A_M4 = A_14 * A_24;
         Vwork = Vwork * A_M4;
         LE24 = patch('faces', FE24, 'vertices' ,Vwork(:,1:3)); 
         set(LE24, 'facec', 1/255*[119, 0, 153]);
         set(LE24, 'EdgeColor','none'); 
         [LXe24,LYe24,LZe24]=ejesF(A_M4,2);
%% CAR
         Vwork = VECR;
         LECR = patch('faces', FECR, 'vertices' ,Vwork(:,1:3)); 
         set(LECR, 'facec', 1/255*[119, 150, 153]);
         set(LECR, 'EdgeColor','none'); 
%% EJES
        setappdata(0, 'A_M1', A_M1);
        setappdata(0, 'A_M2', A_M2);
        setappdata(0, 'A_M3', A_M3);
        setappdata(0, 'A_M4', A_M4);
%% PLOTERS 
        setappdata(0, 'LE21', LE21);
        setappdata(0, 'LE22', LE22);
        setappdata(0, 'LE23', LE23);
        setappdata(0, 'LE24', LE24);
        setappdata(0, 'LCC', LCC); 
        setappdata(0, 'LECR', LECR);
        
    end  %initHome

    %%
    function ejes
        light % add a default light
        daspect([1 1 1]) % Setting the aspect ratio, hace que los círculos se vean como círculos 
        view(135,25)
        xlabel('X'),ylabel('Y'),zlabel('Z');
        
        axis([-80 80 -80 80 -80 80]/3);
        X=[0,80];
        Y=[0,0];
        Z=[0,0];
        line (X,Y,Z,'color','r')
        X=[0,0];
        Y=[0,80];
        Z=[0,0];
        line (X,Y,Z,'color','g')
        X=[0,0];
        Y=[0,0];
        Z=[0,80];
        line (X,Y,Z,'color','b')
        
        plot3([-80,80],[-80,-80],[-80,-80],'k')
        plot3([-80,-80],[-80,80],[-80,-80],'k')
        plot3([-80,-80],[-80,-80],[-80,80],'k')
        plot3([-80,-80],[80,80],[-80,80],'k')
        plot3([-80,80],[-80,-80],[80,80],'k')
        plot3([-80,-80],[-80,80],[80,80],'k')
        grid
    end%ejes

    function[LXe,LYe,LZe]=ejesF(A_01,L)
        bx = A_01(4,1);
        by = A_01(4,2);
        bz = A_01(4,3);
            
       % Toda la transformacion soobre el eje de coordenadas se puede
            % rescatar de la matriz de transformacion A_01
            % EJE X
            VN = [bx by bz]+L*A_01(1,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            LXe = line(EXeX,EXeY,EXeZ, 'color','r');
           % EJE Y
            VN = [bx by bz]+L*A_01(2,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            LYe = line(EXeX,EXeY,EXeZ, 'color','g');
            
            % EJE Z
            VN = [bx by bz]+L*A_01(3,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            LZe = line(EXeX,EXeY,EXeZ, 'color','b');
            
    end%ejesF

    function ejesFA(A_01,L)
        LXe3 = getappdata(0, 'LXe3');
        LYe3 = getappdata(0, 'LYe3');
        LZe3 = getappdata(0, 'LZe3');
        bx = A_01(4,1);
        by = A_01(4,2);
        bz = A_01(4,3);
            
       % Toda la transformacion soobre el eje de coordenadas se puede
            % rescatar de la matriz de transformacion A_01
            % EJE X
            VN = [bx by bz]+L*A_01(1,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            set(LXe3, 'xdata', EXeX, 'ydata', EXeY, 'zdata', EXeZ);
           % EJE Y
            VN = [bx by bz]+L*A_01(2,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            set(LYe3, 'xdata', EXeX, 'ydata', EXeY, 'zdata', EXeZ);
            
            % EJE Z
            VN = [bx by bz]+L*A_01(3,1:3);
            EXeX = [bx VN(1)];
            EXeY = [by VN(2)];
            EXeZ = [bz VN(3)];
            set(LZe3, 'xdata', EXeX, 'ydata', EXeY, 'zdata', EXeZ);
            
    end%ejesFA

    %% Funciones auxiliares
    function Tp = transl(x, y, z)
        t = [x y z];
        Tp = [eye(3)	 zeros(3,1);
            t	1];%Traslación cuando las coordenadas son filas
    end % Funcion de TP

    %% funcion de rotacion Troty
    function T=Troty(phi)%phi se define en grados y radianes
        phi=(pi/180)*phi;
        cp=cos(phi);
        sp=sin(phi);
     % mat transpuesta   
        T=[cp 0 -sp 0;
           0 1 0 0;
           sp 0 cp 0;
           0 0 0 1];
    end %Troty

    %% funcion de rotacion Trotz
    function T=Trotz(phi)%phi se define en grados y radianes
        phi=(pi/180)*phi;
        cp=cos(phi);
        sp=sin(phi);
     % mat transpuesta   
        T=[cp sp 0 0;
           -sp cp 0 0;
           0 0 1 0;
           0 0 0 1];
    end %Troty

    %% funcion de rotacion Trotx
    function T=Trotx(phi)%phi se define en grados y radianes
        phi=(pi/180)*phi;
        cp=cos(phi);
        sp=sin(phi);
     % mat transpuesta   
        T=[1 0 0 0;
           0 cp sp 0;
           0 -sp cp 0;
           0 0 0 1];
    end %Trotx

    %% funcion TDH
    function T = TDH(theta,d,a,alpha)
        theta = theta * pi/180; %conversion a radianes
        alpha = alpha * pi/180; %conversion a radianes
        ct = cos(theta);
        st = sin(theta);
        ca = cos(alpha);
        sa = sin(alpha);
        T = [ct,-ca*st,sa*st,a*ct;
            st,ca*ct,-sa*ct,a*st;
            0,sa,ca,d;
            0,0,0,1];
        T = T';
    end %TDH
    
    %%
    function [q1, q2] = InvKin(xf, yf)
        l1 = getappdata(0, 'L1');
        l2 = getappdata(0, 'L4');
        
        r2 = xf*xf + yf*yf;
        cq = (r2-l1*l1-l2*l2)/(2*l1*l2);
        sq = +sqrt(abs(r2));
        q2 = atan2(sq, cq);
        q2 = q2*180/pi; %conversion a grados
        
        beta = atan2(yf, xf);
        alpha = atan2(sq, cq);
        q1 = beta - alpha;
        q1 = q1*180/pi;%conversion a grados
    end %InvKin

    function [q1, q2, q3] = InvKin3(xf, yf, theta)
        l3 = getappdata(0, 'L1');
        theta = theta*pi/180;
        xff = xf -l3*cos(theta);
        yff = yf -l3*sin(theta);
        [q1, q2] = InvKin(xff, yff);
        theta = theta*180/pi;
        q3 = theta - q1 -q2 - 90;
    end %InvKin
   
    %% Animacion
    function robot_animado(q1fin, q2fin, n)

        qOld = getappdata(0,'qOld');
        l1 = getappdata(0, 'L1');
        l2 = getappdata(0, 'L4');
        offset = getappdata(0, 'offset');

        VE21 = getappdata(0, 'VE21');
        FE21 = getappdata(0, 'FE21');
        VE22 = getappdata(0, 'VE22');
        FE22 = getappdata(0, 'FE22');
        VE23 = getappdata(0, 'VE23');
        FE23 = getappdata(0, 'FE23');
        VE24 = getappdata(0, 'VE24');
        FE24 = getappdata(0, 'FE24');
        VCC = getappdata(0, 'VCC');
        VECR = getappdata(0, 'VECR');
        
        LECR = getappdata(0, 'LECR');
        LE21 = getappdata(0, 'LE21');
        LE22 = getappdata(0, 'LE22');
        LE23 = getappdata(0, 'LE23');
        LE24 = getappdata(0, 'LE24');
        LCC = getappdata(0, 'LCC');
        
        A_M1 = getappdata(0, 'A_M1');
        A_M2 = getappdata(0, 'A_M2');
        A_M3 = getappdata(0, 'A_M3');
        A_M4 = getappdata(0, 'A_M4');
        
        q1Ini = qOld(1);
        II = linspace(q1Ini,q1fin,n);
        for ii=II
            A_01 = transl(0,0,ii);
            VWH = VE21;
            A_M11 = A_M1 * A_01;
            VWH = VWH * A_M11;
            set(LE21, 'vertices', VWH(:,1:3));
            
            VWH = VE21;
            A_M22 = A_M2 * A_01;
            VWH = VWH * A_M22;
            set(LE22, 'vertices', VWH(:,1:3));
            
            VWH = VE21;
            A_M33 = A_M3 * A_01;
            VWH = VWH * A_M33;
            set(LE23, 'vertices', VWH(:,1:3));
            
            VWH = VE21;
            A_M44 = A_M4 * A_01;
            VWH = VWH * A_M44;
            set(LE24, 'vertices', VWH(:,1:3));
            
            VWH = VCC;
            VWH = VWH * A_01;
            set(LCC, 'vertices', VWH(:,1:3));
            
            VWH = VECR;
            VWH = VWH * A_01;
            set(LECR, 'vertices', VWH(:,1:3));
        end
        
        qOld(1) = q1fin;
        setappdata(0,'qOld',qOld);
    end %robot_animado

    function robot_animado2(q1fin, q2fin, n)

        qOld = getappdata(0,'qOld');
        l1 = getappdata(0, 'L1');
        l2 = getappdata(0, 'L4');
        offset = getappdata(0, 'offset');

        VECR = getappdata(0, 'VWCR');
        LECR = getappdata(0, 'LECR');
        
        q1Ini = qOld(1);
        II = linspace(q1Ini,q1fin,n);
        for ii=II
            A_01 = transl(0,-ii,0);
            VWH = v;
            VWH = VWH * A_01;
            set(LECR, 'vertices', VWH(:,1:3));
        end
        qOld(1) = q1fin;
        setappdata(0,'qOld',qOld);
    end %robot_animado

    %% Slider Callback
    function q1_slider_button_press(h,dummy)
        q1fin=get(h, 'value');
        qOld = getappdata(0, 'qOld');
        q2fin = qOld(2);
        robot_animado(q1fin, q2fin ,n)
    end %q1

    function Push_button_press(h,dummy)
        %car_animado();
        %activate = ~activate;
    end%Push_button_press

end%RobotPlanar3