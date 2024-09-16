function [X,Y] = model_code(mode, param, TNF0,IL60,IL100,Temp0,BP0,HR0)
    [tstim1,tstim2,dose,simtime] = input_func(mode);

    k1 = 0; kP = 6.5; k2 = 0;

    % Initial conditions
    y0 = [0; 0; 0; 0; TNF0; 0; IL60; 0; IL100; 0; 0; Temp0; BP0; HR0; 0; 0];
    tstart = 0;
    tfinal = simtime;
    options = odeset('Events',@events,'OutputSel',1,'Refine',4);

    tout = tstart;
    yout = y0.';
    teout = [];
    yeout = [];
    ieout = [];

    while tout<tfinal
        [t,y,te,ye,ie] = ode45(@system,[tstart,tfinal],y0,options);
        nt = length(t);
        if nt<2, break; end
        tout = [tout; t(2:nt)];
        yout = [yout; y(2:nt,:)];
        teout = [teout; te];
        yeout = [yeout; ye];
        ieout = [ieout; ie];
        y0 = y(nt,:);

        if (dose ~= 0 && (t(end)>tstim1 && t(end)<tstim2))
            dose=0.1035;
        elseif (dose ~= 0 && t(end)>tstim2)
            dose=0;
        elseif (k1 ~= 0 && t(end)>47)
            k1=0.01;
        end
        options = odeset(options,'InitialStep',t(nt)-t(nt-4),'MaxStep',t(nt)-t(1));
        tstart = t(nt);
    end

    X = 0:0.01:tout(end);
    Y = interp1(tout,yout,X,'spline');

    function dYdt = system(t,Y)
        % Cytokine scaling (4)
        sTNF = param(1);
        sIL6 = param(2);
        sIL10 = param(3);
        sIL1 = 450;

        % mRNA elimination (4)
        kTNFmRNA = param(4);
        kIL6mRNA = param(5);
        kIL10mRNA= param(6);
        kIL1mRNA = 0.38*0.693;

        % Cytokine elimination (4)
        kTNF = 2*0.693;
        kIL6 = 6*0.693;
        kIL10= 30*0.693;
        kIL1 = 17.85*0.693;

        % Hill equation parameters (18)
        xLPS_TNF = 13.8;
        nLPS_TNF = 1;
        xLPS_IL6 = 0.0006;
        nLPS_IL6 = 0.4;
        xLPS_IL10 = 29.95;
        nLPS_IL10 = 0.64;
        xLPS_IL1 = 0.85;
        nLPS_IL1 = 0.98;
        xIL10_TNF = 0.6;
        nIL10_TNF = 0.9;
        xIL10_IL6 = 15;
        nIL10_IL6 = 0.3;
        xIL10_IL1 = 0.01;
        nIL10_IL1 = 2;
        xIL10_IL10 = 100;
        nIL10_IL10 = 2;
        xIL6_Temp = 3000;
        nIL6_Temp = 2;

        % Heart rate (2)
        kHRBP = 0.5;
        kHRTemp = 15;

        % Differential equation parameters
        LPS_dose = 165;
        kM = 0.33;
        kM1 = 0.1;
        nTemp = 2.7;
        sTempIL1 = 0.01;
        kTemp = 0.9;
        sTempIL6 = 5;
        TempMax = 42;
        sD = 0.0001;
        kD = 0.01;

        tmaxIL6 = 2;
        ntIL6 = 2;
        tmaxIL10 = 2.5;
        ntIL10 = 3;
        tmaxIL1 = 1;
        ntIL1 = 1;


        % Assign state variables
        LPS = Y(1); M = Y(2); M1 = Y(3); mrnaTNF = Y(4); TNF = Y(5);
        mrnaIL6 = Y(6); IL6 = Y(7); mrnaIL10 = Y(8); IL10 = Y(9);
        mrnaIL1 = Y(10); IL1 = Y(11); Temp = Y(12); BP = Y(13); HR = Y(14);
        D = Y(15); dTemp0 = Y(16);

        dP = k1 + dose*LPS_dose - (kP + k2) * LPS;
        dM = LPS - kM * M;
        dM1 = kM * M - kM1 * M1;
        dmrnaTNF = Hup(LPS,xLPS_TNF,nLPS_TNF) * Hdown(IL10,xIL10_TNF,nIL10_TNF) - kTNFmRNA * mrnaTNF;
        dTNF = sTNF * mrnaTNF * M1 - kTNF * (TNF);
        dmrnaIL6 = Hup(LPS,xLPS_IL6,nLPS_IL6) * Hdown(IL10,xIL10_IL6,nIL10_IL6) - kIL6mRNA * mrnaIL6;
        dIL6 = sIL6 * mrnaIL6 * M1 * Hup(t,tmaxIL6,ntIL6) - kIL6 * (IL6);
        dmrnaIL10 = Hup(LPS,xLPS_IL10,nLPS_IL10) * Hdown(IL10,xIL10_IL10,nIL10_IL10) - kIL10mRNA * mrnaIL10;
        dIL10 = sIL10 * mrnaIL10 * Hup(t,tmaxIL10,ntIL10) * M1 - kIL10 * (IL10);
        dmrnaIL1 = Hup(LPS,xLPS_IL1,nLPS_IL1) * Hdown(IL10,xIL10_IL1,nIL10_IL1) - kIL1mRNA * mrnaIL1;
        dIL1 = sIL1 * mrnaIL1 * M1 * Hup(t,tmaxIL1,ntIL1) - kIL1 * (IL1);
        dTemp0 = sTempIL1 * (TempMax - Temp)^nTemp * IL1 + sTempIL6 * Hup(IL6,xIL6_Temp,nIL6_Temp) - kTemp * dTemp0;
        dTemp = kTemp * dTemp0;
        dD = sD * (TNF + IL1 + IL6) - kD * D;
        dBP = -dD;
        dHR = -kHRBP * dBP + kHRTemp * dTemp;

        dYdt = [dP; dM; dM1; dmrnaTNF; dTNF; dmrnaIL6; dIL6; dmrnaIL10; dIL10; dmrnaIL1; dIL1; dTemp; dBP; dHR; dD; dTemp0];
    end

    function Hup = Hup(X,n,h)
        Hup = X^h/(n^h + X^h);
    end

    function Hdown = Hdown(X,n,h)
        Hdown = n^h/(n^h + X^h);
    end

    function [value,isterminal,direction] = events(t,y)
        if (dose==1 && (t>tstim1 && t<tstim2)) || (dose~= 0 && t>tstim2) || (k1~= 0 && t>48)
            value = 0;
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
    end
end
