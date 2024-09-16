function AUC = PSA_model_code(mode, param, TNF0,IL60,IL100,Temp0,BP0,HR0)
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

    AUC = trapz(X,real(Y));

    function dYdt = system(t,Y)
        % Cytokine scaling (4)
        sTNF = param(1);
        sIL6 = param(2);
        sIL10 = param(3);
        sIL1 = param(4);

        % mRNA elimination (4)
        kTNFmRNA = param(5);
        kIL6mRNA = param(6);
        kIL10mRNA = param(7);
        kIL1mRNA = param(8);

        % Cytokine elimination (4)
        kTNF = param(9);
        kIL6 = param(10);
        kIL10 = param(11);
        kIL1 = param(12);

        % Hill equation parameters (18)
        xLPS_TNF = param(13);
        nLPS_TNF = param(14);
        xLPS_IL6 = param(15);
        nLPS_IL6 = param(16);
        xLPS_IL10 = param(17);
        nLPS_IL10 = param(18);
        xLPS_IL1 = param(19);
        nLPS_IL1 = param(20);
        xIL10_TNF = param(21);
        nIL10_TNF = param(22);
        xIL10_IL6 = param(23);
        nIL10_IL6 = param(24);
        xIL10_IL1 = param(25);
        nIL10_IL1 = param(26);
        xIL10_IL10 = param(27);
        nIL10_IL10 = param(28);
        xIL6_Temp = param(29);
        nIL6_Temp = param(30);

        % Heart rate (2)
        kHRBP = param(31);
        kHRTemp = param(32);

        % Differential equation parameters
        LPS_dose = param(33);
        kM = param(34);
        kM1 = param(35);
        nTemp = param(36);
        sTempIL1 = param(37);
        kTemp = param(38);
        sTempIL6 = param(39);
        TempMax = param(40);
        sD = param(41);
        kD = param(42);

        tmaxIL6 = param(43);
        ntIL6 = param(44);
        tmaxIL10 = param(45);
        ntIL10 = param(46);
        tmaxIL1 = param(47);
        ntIL1 = param(48);


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
