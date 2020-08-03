
% Generate coefficients for export profiles fit to a Power Law curve
% From Lauderdale & Cael

b = 0.84; % b-value
zmod = [50 120 220 360 550 790 1080 1420 1810 2250 2740 3280 3870 4510]; % interface depths

z=linspace(50,5200,2576);
Fmar = (z./z(1)).^(-b); % reference power law

[pabsexp,pabsbal,pabsdbl,pabsstr,pabsrat,pabsgam,...
 rabsexp,rabsbal,rabsdbl,rabsstr,rabsrat,rabsgam]=...
    calc_profile_coeffs(zmod,b,'fitabs');

Fabsexp = pabsexp(1)*exp(-z./pabsexp(2));
Fabsbal = pabsbal(1)*(exp(-z./pabsbal(2))+pabsbal(3));
Fabsdbl = pabsdbl(1)*(exp(-z./pabsdbl(2))+pabsdbl(3)*exp(-z./pabsdbl(4)));
Fabsstr = pabsstr(1)*exp(-z.^pabsstr(2));
Fabsrat = pabsrat(1)./(z+pabsrat(2));
Fabsgam = pabsgam(1)*igamma(0,z./pabsgam(2));

coeff_array=zeros(19,5);
coeff_array(1,1  ) = b;
coeff_array(2,1:2) = pabsexp;
coeff_array(3,1:3) = pabsbal;
coeff_array(4,1:4) = pabsdbl;
coeff_array(5,1:2) = pabsstr;
coeff_array(6,1:2) = pabsrat;
coeff_array(7,1:2) = pabsgam;

coeff_array(2,5)   = ser(rabsexp,2);
coeff_array(3,5)   = ser(rabsbal,3);
coeff_array(4,5)   = ser(rabsdbl,4);
coeff_array(5,5)   = ser(rabsstr,2);
coeff_array(6,5)   = ser(rabsrat,2);
coeff_array(7,5)   = ser(rabsgam,2);

[prelexp,prelbal,preldbl,prelstr,prelrat,prelgam,...
 rrelexp,rrelbal,rreldbl,rrelstr,rrelrat,rrelgam]=...
    calc_profile_coeffs(zmod,b,'fitrel');

Frelexp = prelexp(1)*exp(-z./prelexp(2));
Frelbal = prelbal(1)*(exp(-z./prelbal(2))+prelbal(3));
Freldbl = preldbl(1)*(exp(-z./preldbl(2))+preldbl(3)*exp(-z./preldbl(4)));
Frelstr = prelstr(1)*exp(-z.^prelstr(2));
Frelrat = prelrat(1)./(z+prelrat(2));
Frelgam = prelgam(1)*igamma(0,z./prelgam(2));

coeff_array(8 ,1:2) = prelexp;
coeff_array(9 ,1:3) = prelbal;
coeff_array(10,1:4) = preldbl;
coeff_array(11,1:2) = prelstr;
coeff_array(12,1:2) = prelrat;
coeff_array(13,1:2) = prelgam;

coeff_array(8,5)    = ser(rrelexp,2);
coeff_array(9,5)    = ser(rrelbal,3);
coeff_array(10,5)   = ser(rreldbl,4);
coeff_array(11,5)   = ser(rrelstr,2);
coeff_array(12,5)   = ser(rrelrat,2);
coeff_array(13,5)   = ser(rrelgam,2);

[pefdexp,pefdbal,pefddbl,pefdstr,pefdrat,pefdgam,...
 refdexp,refdbal,refddbl,refdstr,refdrat,refdgam]=...
    calc_profile_coeffs(zmod,b,'fitefd');

Fefdexp = pefdexp(1)*exp(-z./pefdexp(2));
Fefdbal = pefdbal(1)*(exp(-z./pefdbal(2))+pefdbal(3));
Fefddbl = pefddbl(1)*(exp(-z./pefddbl(2))+pefddbl(3)*exp(-z./pefddbl(4)));
Fefdstr = pefdstr(1)*exp(-z.^pefdstr(2));
Fefdrat = pefdrat(1)./(z+pefdrat(2));
Fefdgam = pefdgam(1)*igamma(0,z./pefdgam(2));

coeff_array(14,1:2) = pefdexp;
coeff_array(15,1:3) = pefdbal;
coeff_array(16,1:4) = pefddbl;
coeff_array(17,1:2) = pefdstr;
coeff_array(18,1:2) = pefdrat;
coeff_array(19,1:2) = pefdgam;

coeff_array(14,5)   = ser(refdexp,2);
coeff_array(15,5)   = ser(refdbal,3);
coeff_array(16,5)   = ser(refddbl,4);
coeff_array(17,5)   = ser(refdstr,2);
coeff_array(18,5)   = ser(refdrat,2);
coeff_array(19,5)   = ser(refdgam,2);

coeff_array

figure
subplot(231)
plot(Fmar,-z,Frelexp,-z,Frelbal,-z,Freldbl,-z,...
    Frelstr,-z,Frelrat,-z,Frelgam,-z,'LineWidth',2)
legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
title('Minimizing Relative Error')
xlabel('Fraction of export')
ylabel('Depth [m]')
subplot(232)
plot(Fmar,-z,Fabsexp,-z,Fabsbal,-z,Fabsdbl,-z,...
    Fabsstr,-z,Fabsrat,-z,Fabsgam,-z,'LineWidth',2)
legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
title('Minimizing Absolute Error')
xlabel('Fraction of export')
ylabel('Depth [m]')
subplot(233)
plot(Fmar,-z,Fefdexp,-z,Fefdbal,-z,Fefddbl,-z,...
    Fefdstr,-z,Fefdrat,-z,Fefdgam,-z,'LineWidth',2)
legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
title('Fixing e-folding depth')
xlabel('Fraction of export')
ylabel('Depth [m]')

subplot(234)
plot(calc_atten(Fmar,z),-z(1:end-1),...
     calc_atten(Frelexp,z),-z(1:end-1),...
     calc_atten(Frelbal,z),-z(1:end-1),...
     calc_atten(Freldbl,z),-z(1:end-1),...
     calc_atten(Frelstr,z),-z(1:end-1),...
     calc_atten(Frelrat,z),-z(1:end-1),...
     calc_atten(Frelgam,z),-z(1:end-1),...
     'LineWidth',2)
legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')
subplot(235)
plot(calc_atten(Fmar,z),-z(1:end-1),...
     calc_atten(Fabsexp,z),-z(1:end-1),...
     calc_atten(Fabsbal,z),-z(1:end-1),...
     calc_atten(Fabsdbl,z),-z(1:end-1),...
     calc_atten(Fabsstr,z),-z(1:end-1),...
     calc_atten(Fabsrat,z),-z(1:end-1),...
     calc_atten(Fabsgam,z),-z(1:end-1),...
     'LineWidth',2)
 legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')
subplot(236)
plot(calc_atten(Fmar,z),-z(1:end-1),...
     calc_atten(Fefdexp,z),-z(1:end-1),...
     calc_atten(Fefdbal,z),-z(1:end-1),...
     calc_atten(Fefddbl,z),-z(1:end-1),...
     calc_atten(Fefdstr,z),-z(1:end-1),...
     calc_atten(Fefdrat,z),-z(1:end-1),...
     calc_atten(Fefdgam,z),-z(1:end-1),...
     'LineWidth',2)
legend('power','exponential','ballast','double exp',...
    'stretched exp','rational','gamma','Location','southeast')
set(gca,'YLim',[-1000 0]) 
xlabel('Attenuation [m^{-1}]')
ylabel('Depth [m]')
orient landscape
print('export_profiles_matlab.pdf','-dpdf','-bestfit')

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
            % note carefully where functional forms i.e. where parentheses are so we don't make same miscommunication make as before!
            fe = @(p,x)(p(1)*exp(-x./p(2)));
            fb = @(p,x)(p(1)*(exp(-x./p(2))+p(3)));
            fr = @(p,x)(p(1)./(x+p(2)));
            fd = @(p,x)(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4)))));
            fs2 = @(p,x)(p(1)*exp(-x.^p(2))); % TWO PARAMETER STRETCHED EXPONENTIAL
            %fs3 = @(p,x)(p(1)*exp(-p(2).*(x.^p(3)))); % replacing three: can explain over skype if you want
            fg = @(p,x)(p(1)*igamma(0,x./p(2)));

            [pe,re,~,~,~,~] = nlinfit(z,f,fe,[1 1000],'weights',w); % fitting
            [pb,rb,~,~,~,~] = nlinfit(z,f,fb,[1 100 .01],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,f,fr,[100 10],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,f,fd,[1 100 .01 100],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,f,fs3,[10000 10 .1],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,f,fs2,[10 .1],'weights',w);
            [pg,rg,~,~,~,~] = nlinfit(z,f,fg,[.1 1000],'weights',w);

        case('fitrel')

            fe = @(p,x)(log(p(1)*exp(-x./p(2))));
            fb = @(p,x)(log(p(1)*(exp(-x./p(2))+p(3))));
            fr = @(p,x)(log(p(1)./(x+p(2))));
            fd = @(p,x)(log(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4))))));
            fs2= @(p,x)(log(p(1)*exp(-x.^p(2))));
            %fs3 = @(p,x)(log(p(1)*exp(-p(2).*(x.^p(3)))));
            fg = @(p,x)(log(p(1)*igamma(0,x./p(2))));

            [pe,re,~,~,~,~] = nlinfit(z,log(f),fe,[1 1000],'weights',w);
            [pb,rb,~,~,~,~] = nlinfit(z,log(f),fb,[1 100 .01],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,log(f),fr,[100 10],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,log(f),fd,[1 100 .01 100],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,log(f),fs3,[10000 10 .1],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,log(f),fs2,[10 .1],'weights',w);
            [pg,rg,~,~,~,~] = nlinfit(z,log(f),fg,[.1 1000],'weights',w);

        case('fitefd')
            efd = 50*exp(1/b); % add point at efd with lots of weight
            z(end+1) = efd;
            f(end+1) = exp(-1);
            w(end+1) = w(1);

            fe = @(p,x)(p(1)*exp(-x./p(2)));
            fb = @(p,x)(p(1)*(exp(-x./p(2))+p(3)));
            fr = @(p,x)(p(1)./(x+p(2)));
            fd = @(p,x)(p(1)*(exp(-x./p(2))+p(3)*exp(-x./abs(p(4)))));
            fs2= @(p,x)(p(1)*exp(-x.^p(2)));
            %fs3 = @(p,x)(p(1)*exp(-p(2).*(x.^p(3))));
            fg = @(p,x)(p(1)*igamma(0,x./p(2)));

            [pe,re,~,~,~,~] = nlinfit(z,f,fe,[1 1000],'weights',w);
            [pb,rb,~,~,~,~] = nlinfit(z,f,fb,[1 100 .01],'weights',w);
            [pr,rr,~,~,~,~] = nlinfit(z,f,fr,[100 10],'weights',w);
            [pd,rd,~,~,~,~] = nlinfit(z,f,fd,[1 100 .01 100],'weights',w);
           %[ps3,rs,~,~,~,~] = nlinfit(z,f,fs3,[10000 10 .1],'weights',w);
            [ps,rs,~,~,~,~]= nlinfit(z,f,fs2,[10 .1],'weights',w);
            [pg,rg,~,~,~,~] = nlinfit(z,f,fg,[.1 1000],'weights',w);
    end

    % make small correction so F(50m) = 1 exactly
    Fe = pe(1)*exp(-z./pe(2));
    Fb = pb(1)*(exp(-z./pb(2))+pb(3));
    Fr = pr(1)./(z+pr(2));
    Fg = pg(1)*igamma(0,z./pg(2));
    Fs = ps(1)*exp(-z.^ps(2));
    %Fs3 = ps3(1)*exp(-ps3(2).*(z.^ps3(3)));
    Fd = pd(1)*(exp(-z./pd(2))+pd(3)*exp(-z./abs(pd(4))));

    pe(1) = pe(1)./Fe(1);
    pd(1) = pd(1)./Fd(1);
    pr(1) = pr(1)./Fr(1);
    ps(1)= ps(1)./Fs(1);
    pg(1) = pg(1)./Fg(1);
    %ps3(1) = ps3(1)./Fs3(1);
    pb(1) = pb(1)./Fb(1);
     
    % Adjust the residuals to match
    re(1)=f(1)-(pe(1)*exp(-z(1)./pe(2)));
    rb(1)=f(1)-(pb(1)*(exp(-z(1)./pb(2))+pb(3)));
    rd(1)=f(1)-(pr(1)./(z(1)+pr(2)));
    rs(1)=f(1)-(ps(1)*exp(-z(1).^ps(2)));
    rr(1)=f(1)-(pr(1)./(z(1)+pr(2)));
    rg(1)=f(1)-(pg(1)*igamma(0,z(1)./pg(2)));
    
    % assign output arguements
    varargout={pe,pb,pd,ps,pr,pg,... % Coefficient output
               re,rb,rd,rs,rr,rg};   % residuals
end

