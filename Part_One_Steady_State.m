clear
clc
%setting constants
g_K= 36; %combined potassium conductance for all K+ Channels (mSeimens/cm^2)
g_Na = 120; %combined sodium conductance for all Na+ channels (mSeimens/cm^2)
g_L = 0.3; %leakage conductance (mSeimens/cm^2)
E_K = -12; %potasium potential (mV)
E_Na = 115;%sodium potential (mV)
E_L = 10.6; %leakage potential (mV)
Vrest = -70; %Resting membrane potential
Cm = 1.0; %membrane capacitance (microfarad/cm^2)
%V(1)=-70; %membrane potential to be at rest to begin with
V=0; %initializing voltage as 0 for itteration below (the 'for loop')

%setting the step size for the itteration
step = 0.01; %step size set to be 0.1ms

%creating an array for all time values (t) that will be stepped through in
%the test
t = 0:step:100;

%sets no injected current such that the cell is in steady state
I(1:numel(t))=[0];


%Gating variable equations
%where tha a(alpha) represents the probablilitya channel is open and
%b(beta) representing the probability a channel is closed
a_n= 0.01 * ((10-V)/(exp((10-V)/10)-1));
b_n = .125*exp(-V/80);


a_m = .1*((25-V)/(exp((25-V)/10)-1));
b_m = 4*exp(-V/18);


a_h = 0.07*exp(-V/20);
b_h = 1/(exp((30-V)/10)+1);

%initializing values
m(1)= a_m/(a_m + b_m);
n(1)= a_n/(a_n + b_n);
h(1)= a_h/(a_h + b_h);

for i = 1:numel(t)-1
    %alpha and beta values are evaluated at each indivicual time point
    %(i) such and these valuse are used to calculate currents and voltages
    %connected to the voltage gated channels
    a_n(i)= 0.01 * ((10-V(i))/(exp((10-V(i))/10)-1));
    b_n(i) = .125*exp(-V(i)/80);
    a_m(i)= .1*((25-V(i))/(exp((25-V(i))/10)-1));
    b_m(i) = 4*exp(-V(i)/18);
    a_h(i) = 0.07*exp(-V(i)/20);
    b_h(i) = 1/(exp((30-V(i))/10)+1);
    
    %current calculations using the 'm' 'n' and 'h' that tell the
    %probablility a certain channel will be open.  These values were
    %discovered by H&H using curve fitting
    I_Na = (m(i)^3) * g_Na * h(i) * (V(i)-E_Na);
    I_K = (n(i)^4) * g_K * (V(i) - E_K);
    I_L = g_L * (V(i)-E_L);
    %In this case there will be no injected Current (I) therefore there
    %will be no contribution by I(i).  Esentially simplifying the equation
    %to I_ion = I_K - I_Na - I_L;  
    
    %further simplification may be done in this case since you know the
    %cell is in euqilibrium (steady state) the equation may be reduced to the following:      
    %0 = I_K - I_Na - I_L
    
    %the equation is left in its entirety for simplicity of testing other
    %injected current scenarios
    I_ion = I(i) - I_K - I_Na - I_L;
    
    %Euler method for first order differential equations
    % euler method: y(n+1) = y(n)n+h*f(t,y(n))

    V(i+1) = V(i) + step*I_ion/Cm;
    n(i+1) = n(i) + step*(a_n(i) * (1-n(i)) - b_n(i)*n(i));
    m(i+1) = m(i) + step*(a_m(i) * (1-m(i)) - b_m(i)*m(i));
    h(i+1) = h(i) + step*(a_h(i) * (1-h(i)) - b_h(i)*h(i));
    
end
%setting setting initial voltage to be Vrest which is -70mV
V = V-70;

%Voltage plot
plot(t,V)
hold on
legend ({'Membrane Voltage'})
xlabel('Time (ms)')
ylabel('Voltage (mV)')
%axis([0 100 -80 0]) 
title('Membrane Voltage In Steady State')
axis([0 100 -80 0])


%conductance plotting

figure
%plotting the K+ conductance over time
p1 = plot (t,g_K*(n.^4))
hold on
%plotting the Na+ conductance over time
p2 = plot (t, g_Na*(m.^3).*h, 'r')
legend ([p1, p2], 'g_K','g_N_a')
ylabel('Conductance (mS/cm^2)')
xlabel('Time(ms)')
title('Conductance In Steady State')
