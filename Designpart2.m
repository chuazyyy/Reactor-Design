function Designpart2
% Produce a function for integral of heat capacity
    function heatcapacity = ICp(T)

    heatcapacity(1) = 19.037*T+9.146/2*10^-2*T^2-1.217/3*10^-5*T^3-8.033/4*10^-9*T^4;
    heatcapacity(2) = 28.142*T+0.167/2*10^-2*T^2+0.537/3*10^-5*T^3-2.220/4*10^-9*T^4;
    heatcapacity(3) = 29.087*T-0.191/2*10^-2*T^2+0.400/3*10^-5*T^3-0.870/4*10^-9*T^4;
    heatcapacity(4) = 19.874*T+5.021/2*10^-2*T^2+1.268/3*10^-5*T^3-11.00/4*10^-9*T^4;
    heatcapacity(5) = 32.217*T+0.192/2*10^-2*T^2+1.055/3*10^-5*T^3-3.593/4*10^-9*T^4;

end

    CpI=ICp(298); %[J/K mol]
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Dependent Variables: (x(1), x(2), x(3), x(4), x(5))=(m, xi1, xi2, P, T)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define initial flowrates [mol/s]:
    n0total=6300; 
    n0CO = n0total*0.2; 
    n0H2 = n0total*0.75; 
    n0CH3OH = n0total*0.02; 
    n0CH4 = n0total*0.02;
    n0H2O = n0total*0.01; 
    
% Define constants:
    H1=-89980; %[J/mol CH3OH]
    H2=-206296; %[J/mol CH4]
    densitycat = 1400; %[kg/M^3]
    R = 8.314; %[J/K mol]  
    Dp = 0.00704; %[m]
    voidage = 0.4;
    
% Define functions for flowrates [mol/s]:
    nCO = @(x) n0CO - x(2) - x(3); 
    nH2 = @(x) n0H2 - 2.*x(2) - 3.*x(3);
    nCH3OH = @(x) n0CH3OH + x(2);
    nCH4 = @(x) n0CH4 + x(3);
    nH2O = @(x) n0H2O + x(3);
    ntotal = @(x) n0total-2*x(2)-2*x(3);

% Define functions for constants of rate equation:
    K = @(x) 10.^((5332./x(5))-12.831); % [bar^-2]
    A=@(x) -1103.6-308.4.*(x(5)./100)+171.6.*(x(5)./100).^2-14.64.*(x(5)./100).^3;
    B=@(x) 92.65-17.95.*(x(5)./100)-0.5265.*(x(5)./100).^2+0.1745.*(x(5)./100).^3;
    C=@(x) 34.295-5.444.*(x(5)./100)-0.6771.*(x(5)./100).^2+0.1088.*(x(5)./100).^3;
    D=@(x) 483.5-28.21.*(x(5)./100)-22.8.*(x(5)./100).^2+2.438.*(x(5)./100).^3;

% Define denominator for rate equation:
    P=@(x) (A(x)+B(x).*(x(4).*(nCO(x)/ntotal(x)))+C(x).*(x(4).*(nH2(x)/ntotal(x)))+D(x).*(x(4).*(nCH3OH(x)/(ntotal(x))))).^3;

% Rate equations:
    r1=@(x) (x(4).*(nCO(x)./ntotal(x))*(x(4).*(nH2(x)./ntotal(x))).^2-(x(4).*(nCH3OH(x)./ntotal(x)))./K(x))./P(x);
    r2=@(x) (0.5*10^15).*exp(-30000./x(5)).*(x(4).*(nCO(x)./ntotal(x)));

% Specify the molecular weights [kg/mol] of the components in a 5 unit vector
    MW = [0.02801,0.00201,0.03201,0.01604,0.01802];

% G is the superficial mass flux = mdot/A, where A is given as 1 [m^2] 
    A = 1;
    G = sum(MW.*[n0CO,n0H2,n0CH3OH,n0CH4,n0H2O])/A;

% input function for the molecular mass of the gas [kg/mol]
    M = @(x)((nCO(x)./ntotal(x)).*(0.02801))  + (nH2(x)./ntotal(x)).*(0.00201) + (nCH3OH(x)./ntotal(x)).*(0.03204) + (nCH4(x)./ntotal(x)).*(0.01604) + (nH2O(x)./ntotal(x)).*(0.01802);
    
% Define the Cp of all components
    CpCH3OH = @(x) [19.037 9.146*10^-2 -1.217*10^-5 -8.033*10^-9]*[1; x(5); x(5).^2; x(5).^3];
    CpCO = @(x) [28.142 0.167*10^-2 0.537*10^-5 -2.220*10^-9]*[1; x(5); x(5).^2; x(5).^3];
    CpH2 = @(x) [29.087 -0.191*10^-2 0.400*10^-5 -0.870*10^-9]*[1; x(5); x(5).^2; x(5).^3];
    CpCH4 = @(x) [19.874 5.021*10^-2 1.268*10^-5 -11.00*10^-9]*[1; x(5); x(5).^2; x(5).^3];
    CpH2O = @(x) [32.217 0.192*10^-2 1.055*10^-5 -3.593*10^-9]*[1; x(5); x(5).^2; x(5).^3];

% Define the integral of Cp of all components
    CpICH3OH =@(x) [19.037 9.146/2*10^-2 -1.217/3*10^-5 -8.033/4*10^-9]*[x(5); x(5).^2; x(5).^3; x(5).^4];
    CpICO = @(x) [28.142 0.167/2*10^-2 0.537/3*10^-5 -2.220/4*10^-9]*[x(5); x(5).^2; x(5).^3; x(5).^4];
    CpIH2 = @(x) [29.087 -0.191/2*10^-2 0.400/3*10^-5 -0.870/4*10^-9]*[x(5); x(5).^2; x(5).^3; x(5).^4];
    CpICH4 = @(x) [19.874 5.021/2*10^-2 1.268/3*10^-5 -11.00/4*10^-9]*[x(5); x(5).^2; x(5).^3; x(5).^4];
    CpIH2O = @(x) [32.217 0.192/2*10^-2 1.055/3*10^-5 -3.593/4*10^-9]*[x(5); x(5).^2; x(5).^3; x(5).^4];

% Differential equations
    dm= 1400*0.6; 
    dxi1= @(x)1000/3600.*r1(x).*densitycat.*(1-voidage);
    dxi2= @(x)1000/3600.*r2(x).*densitycat.*(1-voidage);
    dp= @(x) 10^-5*(-1.75.*(G^2).*(1-voidage).*R.*x(5))./(Dp.*x(4)*10^5.*M(x).*(voidage^3));
    dT =@(x) (840*1000/3600)*((x(4).*(nCO(x)./ntotal(x))*(x(4).*(nH2(x)/ntotal(x))).^2-(x(4).*(nCH3OH(x)/ntotal(x)))./K(x))./P(x).*(-H1+CpICO(x)+2*CpIH2(x)-CpICH3OH(x)-CpI(2)-2*CpI(3)+CpI(1))+(0.5*10^15).*exp(-30000/x(5)).*(x(4).*((n0CO - x(2) - x(3))/(n0total-2*x(2)-2*x(3)))).*(-H2+CpICO(x)+3*CpIH2(x)-CpICH4(x)-CpIH2O(x)-CpI(2)-3*CpI(3)+CpI(4)+CpI(5)))/((n0CO-x(2)-x(3)).*(CpCO(x))+(n0H2-2.*x(2)-3.*x(3)).*(CpH2(x))+(n0CH3OH+x(2)).*(CpCH3OH(x))+(n0CH4+x(3)).*(CpCH4(x))+(n0H2O+x(3)).*(CpH2O(x)));
    
    
% Solve the ODE system:
    [z,x]=ode45(@(z,x) [dm; dxi1(x); dxi2(x); dp(x); dT(x)], [0 13.8], [0 0 0 225 620]);

% Redefine total amount of individual components based on xi1 and xi2 obtained:
    nCO = @(x) n0CO - x(:,2) - x(:,3);
    nH2 = @(x) n0H2 - 2.*x(:,2) - 3.*x(:,3);
    nCH3OH = @(x) n0CH3OH + x(:,2);
    nCH4 = @(x) n0CH4 + x(:,3);
    nH2O = @(x) n0H2O + x(:,3);
    ntotal = @(x) n0total-2*x(:,2)-2*x(:,3);

% Calculate composition:
    y=[nCO(x)./ntotal(x), nH2(x)./ntotal(x), nCH3OH(x)./ntotal(x), nCH4(x)./ntotal(x), nH2O(x)./ntotal(x)];

% Redefine constants for r1 based on P and T obtained:
    K = @(x) 10.^((5332./x(:,5))-12.831);
    A=@(x) -1103.6-308.4.*(x(:,5)./100)+171.6.*(x(:,5)./100).^2-14.64.*(x(:,5)./100).^3;
    B=@(x) 92.65-17.95.*(x(:,5)./100)-0.5265.*(x(:,5)./100).^2+0.1745.*(x(:,5)./100).^3;
    C=@(x) 34.295-5.444.*(x(:,5)./100)-0.6771.*(x(:,5)./100).^2+0.1088.*(x(:,5)./100).^3;
    D=@(x) 483.5-28.21.*(x(:,5)./100)-22.8.*(x(:,5)./100).^2+2.438.*(x(:,5)./100).^3;
    P= @(x)(A(x)+B(x).*(x(:,4).*(nCO(x)./ntotal(x)))+C(x).*(x(:,4).*(nH2(x)./ntotal(x)))+D(x).*(x(:,4).*(nCH3OH(x)./(ntotal(x))))).^3;

% Calculate r1 based on redefined constants, P, and T:
    rone= (x(:,4).*(nCO(x)./ntotal(x)).*(x(:,4).*(nH2(x)./ntotal(x))).^2-(x(:,4).*(nCH3OH(x)./ntotal(x)))./K(x))./P(x);

% Plot the graphs:
    figure(1);
    plot(z,x(:,2));
    xlabel('Reactor Length');
    ylabel('\xi_1');

    figure(2);
    plot(z,y);
    xlabel('Reactor Length');
    ylabel('Composition');

    figure(3);
    plot(z,x(:,4));
    xlabel('Reactor Length');
    ylabel('Pressure');

    figure(4);
    plot(z,x(:,5));
    xlabel('Reactor Length');
    ylabel('Temperature');

    figure(5);
    plot(x(:,2),x(:,1));
    xlabel('\xi_1');
    ylabel('Catalyst Mass');

    figure(6);
    plot(x(:,2),x(:,3));
    xlabel('\xi_1');
    ylabel('\xi_2');

    figure(7);
    plot(x(:,2),x(:,4));
    xlabel('\xi_1');
    ylabel('Pressure');

    figure(8);
    plot(x(:,2),x(:,5));
    xlabel('\xi_1');
    ylabel('Temperature');

    figure(9);
    plot(x(:,2),rone);
    xlabel('\xi_1');
    ylabel('Rate of Reaction 1');

end