
% Generate coefficients for export profiles fit to a Power Law curve
% From Lauderdale & Cael

b = 0.84; % b-value
zmod = [50 120 220 360 550 790 1080 1420 1810 2250 2740 3280 3870 4510]; % interface depths

z=linspace(50,5200,2576);
Fmar = (z./z(1)).^(-b); % reference power law

[pabsexp,pabsbal,pabsdbl,pabsstr,pabsrat,pabsgam0,pabsgam1,pabsgam2,pabsgam3,pabsgamn,...
 rabsexp,rabsbal,rabsdbl,rabsstr,rabsrat,rabsgam0,rabsgamn]=...
    calc_profile_coeffs(zmod,b,'fitabs');

Fabsexp = pabsexp(1)*exp(-z./pabsexp(2));
Fabsbal = pabsbal(1)*(exp(-z./pabsbal(2))+pabsbal(3));
Fabsdbl = pabsdbl(1)*(exp(-z./pabsdbl(2))+pabsdbl(3)*exp(-z./pabsdbl(4)));
Fabsstr = pabsstr(1)*exp(-z.^pabsstr(2));
Fabsrat = pabsrat(1)./(z+pabsrat(2));
Fabsgam0 = pabsgam0(1)*igamma(0,z./pabsgam0(2));
Fabsgam1 = pabsgam1(1)*igamma(1,z./pabsgam1(2));
Fabsgam2 = pabsgam2(1)*igamma(2,z./pabsgam2(2));
Fabsgam3 = pabsgam3(1)*igamma(3,z./pabsgam3(2));
Fabsgamn = pabsgamn(1)*igamma(pabsgamn(3),z./pabsgamn(2));

coeff_array=zeros(19,5);
coeff_array(1,1  ) = b;
coeff_array(2,1:2) = pabsexp;
coeff_array(3,1:3) = pabsbal;
coeff_array(4,1:4) = pabsdbl;
coeff_array(5,1:2) = pabsstr;
coeff_array(6,1:2) = pabsrat;
coeff_array(7,1:2) = pabsgam0;

coeff_array(2,5)   = ser(rabsexp,2);
coeff_array(3,5)   = ser(rabsbal,3);
coeff_array(4,5)   = ser(rabsdbl,4);
coeff_array(5,5)   = ser(rabsstr,2);
coeff_array(6,5)   = ser(rabsrat,2);
coeff_array(7,5)   = ser(rabsgam0,2);

[prelexp,prelbal,preldbl,prelstr,prelrat,prelgam0,prelgam1,prelgam2,prelgam3,prelgamn,...
 rrelexp,rrelbal,rreldbl,rrelstr,rrelrat,rrelgam0,rrelgamn]=...
    calc_profile_coeffs(zmod,b,'fitrel');

Frelexp = prelexp(1)*exp(-z./prelexp(2));
Frelbal = prelbal(1)*(exp(-z./prelbal(2))+prelbal(3));
Freldbl = preldbl(1)*(exp(-z./preldbl(2))+preldbl(3)*exp(-z./preldbl(4)));
Frelstr = prelstr(1)*exp(-z.^prelstr(2));
Frelrat = prelrat(1)./(z+prelrat(2));
Frelgam0 = prelgam0(1)*igamma(0,z./prelgam0(2));
Frelgam1 = prelgam1(1)*igamma(1,z./prelgam1(2));
Frelgam2 = prelgam2(1)*igamma(2,z./prelgam2(2));
Frelgam3 = prelgam3(1)*igamma(3,z./prelgam3(2));
Frelgamn = prelgamn(1)*igamma(prelgamn(3),z./prelgamn(2));

coeff_array(8 ,1:2) = prelexp;
coeff_array(9 ,1:3) = prelbal;
coeff_array(10,1:4) = preldbl;
coeff_array(11,1:2) = prelstr;
coeff_array(12,1:2) = prelrat;
coeff_array(13,1:2) = prelgam0;

coeff_array(8 ,5)   = ser(rrelexp,2);
coeff_array(9 ,5)   = ser(rrelbal,3);
coeff_array(10,5)   = ser(rreldbl,4);
coeff_array(11,5)   = ser(rrelstr,2);
coeff_array(12,5)   = ser(rrelrat,2);
coeff_array(13,5)   = ser(rrelgam0,2);

[pefdexp,pefdbal,pefddbl,pefdstr,pefdrat,pefdgam0,pefdgam1,pefdgam2,pefdgam3,pefdgamn,...
 refdexp,refdbal,refddbl,refdstr,refdrat,refdgam0,refdgamn]=...
    calc_profile_coeffs(zmod,b,'fitefd');

Fefdexp = pefdexp(1)*exp(-z./pefdexp(2));
Fefdbal = pefdbal(1)*(exp(-z./pefdbal(2))+pefdbal(3));
Fefddbl = pefddbl(1)*(exp(-z./pefddbl(2))+pefddbl(3)*exp(-z./pefddbl(4)));
Fefdstr = pefdstr(1)*exp(-z.^pefdstr(2));
Fefdrat = pefdrat(1)./(z+pefdrat(2));
Fefdgam0 = pefdgam0(1)*igamma(0,z./pefdgam0(2));
Fefdgam1 = pefdgam1(1)*igamma(1,z./pefdgam1(2));
Fefdgam2 = pefdgam2(1)*igamma(2,z./pefdgam2(2));
Fefdgam3 = pefdgam3(1)*igamma(3,z./pefdgam3(2));
Fefdgamn = pefdgamn(1)*igamma(pefdgamn(3),z./pefdgamn(2));

coeff_array(14,1:2) = pefdexp;
coeff_array(15,1:3) = pefdbal;
coeff_array(16,1:4) = pefddbl;
coeff_array(17,1:2) = pefdstr;
coeff_array(18,1:2) = pefdrat;
coeff_array(19,1:2) = pefdgam0;

coeff_array(14,5)   = ser(refdexp,2);
coeff_array(15,5)   = ser(refdbal,3);
coeff_array(16,5)   = ser(refddbl,4);
coeff_array(17,5)   = ser(refdstr,2);
coeff_array(18,5)   = ser(refdrat,2);
coeff_array(19,5)   = ser(refdgam0,2);

coeff_array

%%
figure
subplot(231); hold on
plot(Frelgamn,-z,'color','#cccccc','LineWidth',2)
plot(Frelgam0,-z,'color','#999999','LineWidth',2)
plot(Frelgam1,-z,'color','#808080','LineWidth',2)
plot(Frelgam2,-z,'color','#4d4d4d','LineWidth',2)
plot(Frelgam3,-z,'color','#000000','LineWidth',2)
plot(Fmar   ,-z,'--','color','#1f67b4','LineWidth',2)
plot(Frelexp,-z,'-.','color','#ff6f0e','LineWidth',2)
plot(Frelbal,-z,'color','#2ca02c','LineWidth',2)
plot(Freldbl,-z,'color','#9467bd','LineWidth',2)
plot(Frelstr,-z,'color','#8c564b','LineWidth',2)
plot(Frelrat,-z,'color','#d62728','LineWidth',2)
plot(Frelgam0,-z,':','color','#e377c2','LineWidth',2)
%legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
%    'power','exponential','ballast','double exp',...
%    'stretched exp','rational','gamma',...
%    'Location','southeast')
set(gca,'XLim',[0 1],'YLim',[-1000 0])
box on
title('Minimizing Relative Error')
xlabel('Fraction of export')
ylabel('Depth [m]')
subplot(232); hold on
plot(Fabsgamn,-z,'color','#cccccc','LineWidth',2)
plot(Fabsgam0,-z,'color','#999999','LineWidth',2)
plot(Fabsgam1,-z,'color','#808080','LineWidth',2)
plot(Fabsgam2,-z,'color','#4d4d4d','LineWidth',2)
plot(Fabsgam3,-z,'color','#000000','LineWidth',2)
plot(Fmar   ,-z,'--','color','#1f67b4','LineWidth',2)
plot(Fabsexp,-z,'-.','color','#ff6f0e','LineWidth',2)
plot(Fabsbal,-z,'color','#2ca02c','LineWidth',2)
plot(Fabsdbl,-z,'color','#9467bd','LineWidth',2)
plot(Fabsstr,-z,'color','#8c564b','LineWidth',2)
plot(Fabsrat,-z,'color','#d62728','LineWidth',2)
plot(Fabsgam0,-z,':','color','#e377c2','LineWidth',2)
legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
    'power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma',...
    'Location','southeast')
set(gca,'XLim',[0 1],'YLim',[-1000 0]) 
box on
title('Minimizing Absolute Error')
xlabel('Fraction of export')
ylabel('Depth [m]')
subplot(233); hold on
plot(Fefdgamn,-z,'color','#cccccc','LineWidth',2)
plot(Fefdgam0,-z,'color','#999999','LineWidth',2)
plot(Fefdgam1,-z,'color','#808080','LineWidth',2)
plot(Fefdgam2,-z,'color','#4d4d4d','LineWidth',2)
plot(Fefdgam3,-z,'color','#000000','LineWidth',2)
plot(Fmar   ,-z,'--','color','#1f67b4','LineWidth',2)
plot(Fefdexp,-z,'-.','color','#ff6f0e','LineWidth',2)
plot(Fefdbal,-z,'color','#2ca02c','LineWidth',2)
plot(Fefddbl,-z,'color','#9467bd','LineWidth',2)
plot(Fefdstr,-z,'color','#8c564b','LineWidth',2)
plot(Fefdrat,-z,'color','#d62728','LineWidth',2)
plot(Fefdgam0,-z,':','color','#e377c2','LineWidth',2)
legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
    'power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma',...
    'Location','southeast')
set(gca,'XLim',[0 1],'YLim',[-1000 0]) 
box on
title('Fixing e-folding depth')
xlabel('Fraction of export')
ylabel('Depth [m]')

subplot(234); hold on
plot(calc_atten(Frelgamn,z),-z(1:end-1),'color','#cccccc','LineWidth',2)
plot(calc_atten(Frelgam0,z),-z(1:end-1),'color','#999999','LineWidth',2)
plot(calc_atten(Frelgam1,z),-z(1:end-1),'color','#808080','LineWidth',2)
plot(calc_atten(Frelgam2,z),-z(1:end-1),'color','#4d4d4d','LineWidth',2)
plot(calc_atten(Frelgam3,z),-z(1:end-1),'color','#000000','LineWidth',2)
plot(calc_atten(Fmar   ,z),-z(1:end-1),'--','color','#1f67b4','LineWidth',2)
plot(calc_atten(Frelexp,z),-z(1:end-1),'-.','color','#ff6f0e','LineWidth',2)
plot(calc_atten(Frelbal,z),-z(1:end-1),'color','#2ca02c','LineWidth',2)
plot(calc_atten(Freldbl,z),-z(1:end-1),'color','#9467bd','LineWidth',2)
plot(calc_atten(Frelstr,z),-z(1:end-1),'color','#8c564b','LineWidth',2)
plot(calc_atten(Frelrat,z),-z(1:end-1),'color','#d62728','LineWidth',2)
plot(calc_atten(Frelgam0,z),-z(1:end-1),':','color','#e377c2','LineWidth',2)
legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
    'power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma',...
    'Location','southeast')
set(gca,'XLim',[0 0.02],'YLim',[-1000 0]) 
box on
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')

subplot(235); hold on
plot(calc_atten(Fabsgamn,z),-z(1:end-1),'color','#cccccc','LineWidth',2)
plot(calc_atten(Fabsgam0,z),-z(1:end-1),'color','#999999','LineWidth',2)
plot(calc_atten(Fabsgam1,z),-z(1:end-1),'color','#808080','LineWidth',2)
plot(calc_atten(Fabsgam2,z),-z(1:end-1),'color','#4d4d4d','LineWidth',2)
plot(calc_atten(Fabsgam3,z),-z(1:end-1),'color','#000000','LineWidth',2)
plot(calc_atten(Fmar   ,z),-z(1:end-1),'--','color','#1f67b4','LineWidth',2)
plot(calc_atten(Fabsexp,z),-z(1:end-1),'-.','color','#ff6f0e','LineWidth',2)
plot(calc_atten(Fabsbal,z),-z(1:end-1),'color','#2ca02c','LineWidth',2)
plot(calc_atten(Fabsdbl,z),-z(1:end-1),'color','#9467bd','LineWidth',2)
plot(calc_atten(Fabsstr,z),-z(1:end-1),'color','#8c564b','LineWidth',2)
plot(calc_atten(Fabsrat,z),-z(1:end-1),'color','#d62728','LineWidth',2)
plot(calc_atten(Fabsgam0,z),-z(1:end-1),':','color','#e377c2','LineWidth',2)
%legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
%    'power','exponential','ballast','double exp',...
%    'stretched exp','rational','gamma',...
%    'Location','southeast')
set(gca,'XLim',[0 0.02],'YLim',[-1000 0]) 
box on
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')

subplot(236); hold on
plot(calc_atten(Fefdgamn,z),-z(1:end-1),'color','#cccccc','LineWidth',2)
plot(calc_atten(Fefdgam0,z),-z(1:end-1),'color','#999999','LineWidth',2)
plot(calc_atten(Fefdgam1,z),-z(1:end-1),'color','#808080','LineWidth',2)
plot(calc_atten(Fefdgam2,z),-z(1:end-1),'color','#4d4d4d','LineWidth',2)
plot(calc_atten(Fefdgam3,z),-z(1:end-1),'color','#000000','LineWidth',2)
plot(calc_atten(Fmar   ,z),-z(1:end-1),'--','color','#1f67b4','LineWidth',2)
plot(calc_atten(Fefdexp,z),-z(1:end-1),'-.','color','#ff6f0e','LineWidth',2)
plot(calc_atten(Fefdbal,z),-z(1:end-1),'color','#2ca02c','LineWidth',2)
plot(calc_atten(Fefddbl,z),-z(1:end-1),'color','#9467bd','LineWidth',2)
plot(calc_atten(Fefdstr,z),-z(1:end-1),'color','#8c564b','LineWidth',2)
plot(calc_atten(Fefdrat,z),-z(1:end-1),'color','#d62728','LineWidth',2)
plot(calc_atten(Fefdgam0,z),-z(1:end-1),':','color','#e377c2','LineWidth',2)
%legend('gamma(a,x)','gamma(0,x)','gamma(1,x)','gamma(2,x)','gamma(3,x)',...
%    'power','exponential','ballast','double exp',...
%    'stretched exp','rational','gamma',...
%    'Location','southeast')
set(gca,'XLim',[0 0.02],'YLim',[-1000 0]) 
box on
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')
set(gcf, 'PaperType', 'A4');
orient landscape
print('export_profiles_matlab.pdf','-dpdf','-bestfit')
%%
% write out to a file
writematrix(coeff_array','export_profile_coefficients.csv')

function a = calc_atten(prof,depth)
    a = (1./prof(1:end-1)).*(diff(prof)./diff(-depth));
end

function s=ser(residuals,nparms)
    % Standard Error of Regression (S)
    s = sqrt(sum(residuals.*residuals)/(length(residuals) - nparms));
end 

function varargout = calc_profile_coeffs(z,varargin)
    
    if nargin==1
        b=0.9;              % default martin exponent
        optfitstr='fitabs';   % default fit (minimizing absolute error)   
    elseif nargin==2
        if ischar(varargin{1})
            b=0.9;
            optfitstr=varargin{1};
        elseif strcmpi(class(varargin{1}),'double')
            b      =varargin{1};
            optfitstr='fitabs';
        else
            error('Could not parse input arguments, maybe you supplied a cell array?')
        end
    elseif nargin==3
        b      =varargin{1};
        optfitstr=varargin{2};
    end

    if ~any(strcmpi({'fitabs','fitrel','fitefd'},optfitstr))
            error('Fit type not recognised, needs to be one of "FITABS", "FITREL", or "FITEFD"')
    end
    
    f = (z./z(1)).^(-b); % reference power law
    w = ones(size(z)); w(1) = 1000; % force through z(1) via weighting it heavily

    switch(lower(optfitstr))
        case ('fitabs')
            % Fit the reference curve
            fe = @(p,x)(p(1)*exp(-x./p(2)));
            fb = @(p,x)(p(1)*(exp(-x./p(2))+p(3)));
            fd = @(p,x)(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4)))));
            fs2 = @(p,x)(p(1)*exp(-x.^p(2))); % TWO PARAMETER STRETCHED EXPONENTIAL
            %fs3 = @(p,x)(p(1)*exp(-p(2).*(x.^p(3)))); % replacing three: can explain over skype if you want
            fr = @(p,x)(p(1)./(x+p(2)));
            fg0 = @(p,x)(p(1)*igamma(0,x./p(2)));
            fg1 = @(p,x)(p(1)*igamma(1,x./p(2)));
            fg2 = @(p,x)(p(1)*igamma(2,x./p(2)));
            fg3 = @(p,x)(p(1)*igamma(3,x./p(2)));
            fgx = @(p,x)(p(1)*igamma(p(3),x./p(2)));
            
            [pe,re,~,~,~,~] = nlinfit(z,f,fe,[1 1000],'weights',w); % fitting
            [pb,rb,~,~,~,~] = nlinfit(z,f,fb,[1 100 .01],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,f,fd,[1 100 .01 100],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,f,fs2,[10 .1],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,f,fs3,[10000 10 .1],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,f,fr,[100 10],'weights',w);
            [pg0,rg0,~,~,~,~] = nlinfit(z,f,fg0,[.1 1000],struct('MaxIter',100),'weights',w);
            [pg1,rg1,~,~,~,~] = nlinfit(z,f,fg1,[1 1000],struct('MaxIter',500),'weights',w);
            [pg2,rg2,~,~,~,~] = nlinfit(z,f,fg2,[1 1000],struct('MaxIter',500),'weights',w);
            [pg3,rg3,~,~,~,~] = nlinfit(z,f,fg3,[1 1000],struct('MaxIter',500),'weights',w);
            [pgx,rgx,~,~,~,~] = nlinfit(z,f,fgx,[.1 1000 1],struct('MaxIter',500),'weights',w);

        case('fitrel')

            fe = @(p,x)(log(p(1)*exp(-x./p(2))));
            fb = @(p,x)(log(p(1)*(exp(-x./p(2))+p(3))));
            fd = @(p,x)(log(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4))))));
            fs2= @(p,x)(log(p(1)*exp(-x.^p(2))));
            %fs3 = @(p,x)(log(p(1)*exp(-p(2).*(x.^p(3)))));
            fr = @(p,x)(log(p(1)./(x+p(2))));
            fg0 = @(p,x)(log(p(1)*igamma(0,x./p(2))));
            fg1 = @(p,x)(log(p(1)*igamma(1,x./p(2))));
            fg2 = @(p,x)(log(p(1)*igamma(2,x./p(2))));
            fg3 = @(p,x)(log(p(1)*igamma(3,x./p(2))));
            fgx = @(p,x)(log(p(1)*igamma(p(3),x./p(2))));

            [pe,re,~,~,~,~] = nlinfit(z,log(f),fe,[1 1000],'weights',w);
            [pb,rb,~,~,~,~] = nlinfit(z,log(f),fb,[1 100 .01],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,log(f),fd,[1 100 .01 100],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,log(f),fs2,[10 .1],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,log(f),fs3,[10000 10 .1],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,log(f),fr,[100 10],'weights',w);
            [pg0,rg0,~,~,~,~] = nlinfit(z,log(f),fg0,[.1 1000],struct('MaxIter',100),'weights',w);
            [pg1,rg1,~,~,~,~] = nlinfit(z,log(f),fg1,[1 1000],struct('MaxIter',500),'weights',w);
            [pg2,rg2,~,~,~,~] = nlinfit(z,log(f),fg2,[1 1000],struct('MaxIter',500),'weights',w);
            [pg3,rg3,~,~,~,~] = nlinfit(z,log(f),fg3,[1 1000],struct('MaxIter',500),'weights',w);
            [pgx,rgx,~,~,~,~] = nlinfit(z,log(f),fgx,[.1 1000 1],struct('MaxIter',500),'weights',w);

        case('fitefd')
            efd = 50*exp(1/b); % add point at efd with lots of weight
            z(end+1) = efd;
            f(end+1) = exp(-1);
            w(end+1) = w(1);

            fe = @(p,x)(p(1)*exp(-x./p(2)));
            fb = @(p,x)(p(1)*(exp(-x./p(2))+p(3)));
            fd = @(p,x)(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4)))));
            fs2= @(p,x)(p(1)*exp(-x.^p(2)));
            %fs3 = @(p,x)(p(1)*exp(-p(2).*(x.^p(3))));
            fr = @(p,x)(p(1)./(x+p(2)));
            fg0 = @(p,x)(p(1)*igamma(0,x./p(2)));
            fg1 = @(p,x)(p(1)*igamma(1,x./p(2)));
            fg2 = @(p,x)(p(1)*igamma(2,x./p(2)));
            fg3 = @(p,x)(p(1)*igamma(3,x./p(2)));
            fgx = @(p,x)(p(1)*igamma(p(3),x./p(2)));

            [pe,re,~,~,~,~] = nlinfit(z,f,fe,[1 1000],'weights',w);
            [pb,rb,~,~,~,~] = nlinfit(z,f,fb,[1 100 .01],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,f,fd,[1 100 .01 100],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,f,fs2,[10 .1],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,f,fs3,[10000 10 .1],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,f,fr,[100 10],'weights',w);
            [pg0,rg0,~,~,~,~] = nlinfit(z,f,fg0,[.1 1000],struct('MaxIter',100),'weights',w);
            [pg1,rg1,~,~,~,~] = nlinfit(z,f,fg1,[1 1000],struct('MaxIter',500),'weights',w);
            [pg2,rg2,~,~,~,~] = nlinfit(z,f,fg2,[10 1000],struct('MaxIter',500),'weights',w);
            [pg3,rg3,~,~,~,~] = nlinfit(z,f,fg3,[1 1000],struct('MaxIter',500),'weights',w);
            [pgx,rgx,~,~,~,~] = nlinfit(z,f,fgx,[.1 1000 1],struct('MaxIter',500),'weights',w);
    end

    % make small correction so F(50m) = 1 exactly
    Fe = pe(1)*exp(-z./pe(2));
    Fb = pb(1)*(exp(-z./pb(2))+pb(3));
    Fd = pd(1)*(exp(-z./pd(2))+pd(3)*exp(-z./abs(pd(4))));
    Fs = ps(1)*exp(-z.^ps(2));
    %Fs3 = ps3(1)*exp(-ps3(2).*(z.^ps3(3)));
    Fr = pr(1)./(z+pr(2));
    Fg0 = pg0(1)*igamma(0,z./pg0(2));
    Fg1 = pg1(1)*igamma(1,z./pg1(2));
    Fg2 = pg2(1)*igamma(2,z./pg2(2));
    Fg3 = pg3(1)*igamma(3,z./pg3(2));
    Fgx = pgx(1)*igamma(pgx(3),z./pgx(2));

    pe(1) = pe(1)./Fe(1);
    pb(1) = pb(1)./Fb(1);    
    pd(1) = pd(1)./Fd(1);
    ps(1)= ps(1)./Fs(1);
    %ps3(1) = ps3(1)./Fs3(1);
    pr(1) = pr(1)./Fr(1);
    pg0(1) = pg0(1)./Fg0(1);
    pg1(1) = pg1(1)./Fg1(1);
    pg2(1) = pg2(1)./Fg2(1);
    pg3(1) = pg3(1)./Fg3(1);
    pgx(1) = pgx(1)./Fgx(1);
     
    % Adjust the residuals to match
    re(1)=f(1)-(pe(1)*exp(-z(1)./pe(2)));
    rb(1)=f(1)-(pb(1)*(exp(-z(1)./pb(2))+pb(3)));
    rd(1)=f(1)-(pr(1)./(z(1)+pr(2)));
    rs(1)=f(1)-(ps(1)*exp(-z(1).^ps(2)));
    rr(1)=f(1)-(pr(1)./(z(1)+pr(2)));
    rg0(1)=f(1)-(pg0(1)*igamma(0,z(1)./pg0(2)));
    rg1(1)=f(1)-(pg1(1)*igamma(1,z(1)./pg1(2)));
    rg2(1)=f(1)-(pg2(1)*igamma(2,z(1)./pg2(2)));
    rg3(1)=f(1)-(pg3(1)*igamma(3,z(1)./pg3(2)));
    rgx(1)=f(1)-(pgx(1)*igamma(pgx(3),z(1)./pgx(2)));
    
    % assign output arguements
    varargout={pe,pb,pd,ps,pr,pg0,pg1,pg2,pg3,pgx,... % Coefficient output
               re,rb,rd,rs,rr,rg0,rgx};   % residuals
end

