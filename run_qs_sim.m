function out = run_qs_sim(ov)

% RUN_QS_SIM  Batch run of the calibrated QS-repressilator PDE. Returns the
% period and amplitude of the edge CPK/N observable.
%
%   out = run_qs_sim()                      % nominal calibrated run
%   out = run_qs_sim(struct('D_G',1.9))     % override one or more parameters
%
% OUTPUT struct:
%   out.period_pk   period (h), peak detection      (primary, matches calibration)
%   out.period_acf  period (h), autocorrelation     (independent cross-check)
%   out.amp         normalized peak-to-valley amplitude
%   out.n_peaks     number of detected peaks         (regime indicator)
%   out.t, out.CN   time (h) and CPK/N edge trace    (for optional plotting)
%   out.blewup      true if the run went non-finite

if nargin < 1 || isempty(ov), ov = struct(); end

% ---------- QS PARAMETERS ----------
P.n=10; P.tau=0.1; P.mu=0.3; P.k=100;
P.K_N=2000; P.D_0=10000; P.eps_border=0.01;
P.D_A=0; P.D_G=1; P.D_R=0; P.D_R2=0; P.D_C=0;
P.alpha_A=1; P.alpha_G=1; P.alpha_R=1.1; P.alpha_R2=1.8; P.alpha_C=3;
P.delta_A=1; P.delta_G=1; P.delta_R=0.1; P.delta_R2=0.5; P.delta_C=0.25;
P.K_R2A=0.4; P.K_RR=0.5; P.K_RC=1; P.K_RR2=0.5; P.K_R2C=1;
P.tau_cal=0.70597;

% ---------- APPLY OVERRIDES ----------
fn = fieldnames(ov);
for i=1:numel(fn)
    if ~isfield(P, fn{i})
        error('run_qs_sim:unknownParam','Unknown parameter: %s', fn{i});
    end
    P.(fn{i}) = ov.(fn{i});
end

% unpack
n=P.n; tau=P.tau; mu=P.mu; k=P.k; K_N=P.K_N; D_0=P.D_0; eps_border=P.eps_border;
D_A=P.D_A; D_G=P.D_G; D_R=P.D_R; D_R2=P.D_R2; D_C=P.D_C;
alpha_A=P.alpha_A; alpha_G=P.alpha_G; alpha_R=P.alpha_R; alpha_R2=P.alpha_R2; alpha_C=P.alpha_C;
delta_A=P.delta_A; delta_G=P.delta_G; delta_R=P.delta_R; delta_R2=P.delta_R2; delta_C=P.delta_C;
K_R2A=P.K_R2A; K_RR=P.K_RR; K_RC=P.K_RC; K_RR2=P.K_RR2; K_R2C=P.K_R2C;
tau_cal=P.tau_cal;

% ---------- RESOLUTION / HORIZON (resolution) ----------
dx_base=1; dt_base=0.1; dx_ref=0.125; T_ref=49000;
dx=0.25; s=dx/dx_base; dt=s^2*dt_base;
tFinal_ref = T_ref*((dx_ref/dx_base)^2*dt_base);
T = max(1, round(tFinal_ref/dt));
USE_LAPLACIAN_9 = true;

% ---------- DOMAIN ----------
Lx=60; Ly=60;
Nx=max(4,round(Lx/dx)); Ny=max(4,round(Ly/dx));

% ---------- INITIAL CONDITIONS ----------
[X,Y]=meshgrid(1:Nx,1:Ny);
x0=Nx/2; y0=Ny/2; cx=round(x0); cy=round(y0);
sigma_phys=3; sigma=sigma_phys/dx;
N = K_N*0.01*exp(-((X-x0).^2+(Y-y0).^2)/(2*sigma^2));
A=zeros(Ny,Nx); G=A; R=A; R2=A; C=A; A(cy,cx)=0.1;
R_cells = hypot(X-cx,Y-cy)*dx;   % length units

if USE_LAPLACIAN_9
    lap = @(U) (0.25*(circshift(U,[1,1])+circshift(U,[-1,-1])+circshift(U,[1,-1])+circshift(U,[-1,1])) + ...
        0.5*(circshift(U,[1,0])+circshift(U,[-1,0])+circshift(U,[0,1])+circshift(U,[0,-1])) - 3*U)/dx^2;
else
    lap = @(U) (circshift(U,[1,0])+circshift(U,[-1,0])+circshift(U,[0,1])+circshift(U,[0,-1]) - 4*U)/dx^2;
end

% ---------- OBSERVABLE STORAGE ----------
N_occ = 0.05*K_N; ring_width = 0.5;
time_hours = (0:T-1)'*dt*tau_cal;
CN_edge = zeros(T,1);

% ---------- TIME LOOP ----------
blewup=false; tiLast=T;
for ti=1:T
    Hgate = 1./(1+exp(-(N-2*K_N)/eps_border));
    D = D_0*Hgate;
    dNdx=(circshift(N,[0,-1])-circshift(N,[0,1]))/(2*dx);
    dNdy=(circshift(N,[-1,0])-circshift(N,[1,0]))/(2*dx);
    dDdx=(circshift(D,[0,-1])-circshift(D,[0,1]))/(2*dx);
    dDdy=(circshift(D,[-1,0])-circshift(D,[1,0]))/(2*dx);
    div_grad = D.*lap(N) + dDdx.*dNdx + dDdy.*dNdy;
    growth = mu*N.*(1-N/K_N);
    N = N + dt*(div_grad+growth);
    N = max(N,0); N = min(N,1.1*K_N);

    lower_threshold=tau/10; upper_threshold=7.5*tau;
    H_G = 1./(1+exp(k*(G-lower_threshold))) + 1./(1+exp(k*(upper_threshold-G)));

    A  = A  + dt*(D_A *lap(A)  + alpha_A *N./(1+(R2/K_R2A).^n)                  - delta_A *A);
    G  = G  + dt*(D_G *lap(G)  + alpha_G *A                                     - delta_G *G);
    R  = R  + dt*(D_R *lap(R)  + alpha_R *N./(1+(R/K_RR ).^n.*H_G)              - delta_R *R);
    R2 = R2 + dt*(D_R2*lap(R2) + alpha_R2*N./(1+(R/K_RR2).^n.*H_G)             - delta_R2*R2);
    C  = C  + dt*(D_C *lap(C)  + alpha_C *N./(1+(R2/K_R2C).^n+(R/K_RC).^n.*H_G) - delta_C *C);

    if any(~isfinite(C(:))) || any(~isfinite(N(:)))
        blewup=true; tiLast=ti; break;
    end

    occ = N > N_occ;
    if any(occ(:))
        front = max(R_cells(occ));
        ring_mask = (R_cells<=front)&(R_cells>=front-ring_width)&occ;
        CN_edge(ti) = mean(C(ring_mask)./N(ring_mask));
    else
        CN_edge(ti) = C(cy,cx)/max(N(cy,cx),N_occ);
    end
end

% ---------- BLOW-UP GUARD ----------
if blewup
    out = struct('period_pk',NaN,'period_acf',NaN,'amp',NaN,'n_peaks',0, ...
                 't',time_hours(1:tiLast),'CN',CN_edge(1:tiLast),'blewup',true);
    return;
end

% ---------- MEASURE PERIOD + AMPLITUDE (matches calibration method) ----------
sig = CN_edge; tt = time_hours; dt_h = tt(2)-tt(1);
sig_n = (sig-min(sig))/(max(sig)-min(sig));
rng = max(sig_n)-min(sig_n);
mind = round(8/dt_h);
[~,pk] = findpeaks( sig_n,'MinPeakProminence',0.2*rng,'MinPeakDistance',mind);
[~,tr] = findpeaks(-sig_n,'MinPeakProminence',0.2*rng,'MinPeakDistance',mind);
if numel(pk)>=2, period_pk = mean(diff(tt(pk))); else, period_pk = NaN; end
ncyc = min(numel(pk),numel(tr));
if ncyc>=1, amp = mean(sig_n(pk(1:ncyc))-sig_n(tr(1:ncyc))); else, amp = NaN; end

% ---------- AUTOCORRELATION CROSS-CHECK ----------
xx = sig - mean(sig);
[ac,lg] = xcorr(xx,'coeff');
keep = lg>=0; ac=ac(keep); lg=lg(keep);
[~,acpk] = findpeaks(ac,'MinPeakProminence',0.05,'MinPeakDistance',mind);
if ~isempty(acpk), period_acf = lg(acpk(1))*dt_h; else, period_acf = NaN; end

out = struct('period_pk',period_pk,'period_acf',period_acf,'amp',amp, ...
             'n_peaks',numel(pk),'t',tt,'CN',sig,'blewup',false);
end
