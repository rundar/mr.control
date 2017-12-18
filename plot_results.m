function plot_results(RF, Gs, dt, d)
%     written by Christoph Aigner, Dec 18, 2017
%############## Copyright: ##################################
%  Copyright (C) 2017  Armin Rund, armin.rund@uni-graz.at
%                Christoph Aigner, christoph.aigner@tugraz.at
% Implementation of the time optimal control method of paper
% "Simultaneous Multislice Refocusing via Time Optimal Control"
% by A. Rund, C. Aigner, K. Kunisch, and R. Stollberger.
%%###########################################################

%Bloch foreward simulation (spin domain)
[~,b] = blochsd(RF, Gs, dt, d.xdis, d.gamma);

%compute variables needed for the plot
T = length(RF)*dt; %overall duration in ms
tdis=0:dt:T-dt;    %time vector in ms
Gs_slew=diff(Gs)/dt; %slew rate of Gs
outslice=double(d.outslice); outslice(d.outslice==0)=nan; %outslice
inslice=double(d.inslice); inslice(d.inslice==0)=nan;     %inslice

%plot RF, Gs, Gs_slew and the refocusing profile abs(b.^2)
figure
subplot(2,2,1); hold on; grid on;
if isfield(d,'RF_max'),RF_max=d.RF_max; else, RF_max=d.u_box; end;
plot(tdis, -RF_max*ones(size(RF,1),1),'k'); 
plot(tdis, RF_max*ones(size(RF,1),1),'k');
plot(tdis, RF); 
axis([0 T -RF_max*1.2 RF_max*1.2]);
xlabel('time in ms');
ylabel('RF in mT');

subplot(2,2,2); hold on; grid on;
if isfield(d,'Gs_max'),Gs_max=d.Gs_max; else, Gs_max=d.w_box; end;
plot(tdis, -Gs_max*ones(size(Gs,1),1),'k'); 
plot(tdis, Gs_max*ones(size(Gs,1),1),'k');
plot(tdis, Gs);
axis([0 T -Gs_max*1.2 Gs_max*1.2]);
xlabel('time in ms');
ylabel('Gs in mT/m');

subplot(2,2,3); hold on; grid on;
if isfield(d,'Gs_slew_max'),Gs_slew_max=d.Gs_slew_max; else, Gs_slew_max=d.slew_box; end;
plot(tdis(1:end-1), -Gs_slew_max*ones(size(Gs_slew,1),1),'k');
plot(tdis(1:end-1), Gs_slew_max*ones(size(Gs_slew,1),1),'k');
plot(tdis(1:end-1), Gs_slew);
axis([0 T -Gs_slew_max*1.2 Gs_slew_max*1.2]);
xlabel('time in ms');
ylabel('slew rate in T/m/s');

subplot(2,2,4);
plot(d.xdis*1000,inslice,'k'); hold all; grid on;
plot(d.xdis*1000,inslice-d.maxErr_in,'k'); hold all; grid on;
plot(d.xdis*1000,1-outslice,'k'); hold all; grid on;
plot(d.xdis*1000,1-outslice+d.maxErr_out,'k'); hold all; grid on;
plot(d.xdis*1000,abs(b.^2)); hold all; grid on;
axis([min(d.xdis*1000),max(d.xdis*1000),-0.1, 1.1]);
xlabel('distance in mm');
ylabel('refocusing profile')