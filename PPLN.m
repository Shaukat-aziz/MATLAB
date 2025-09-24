% PPLN_designer_interactive_v2.m
% Interactive PPLN designer v2
% - JSI + phase-matching visualization
% - Monte-Carlo sampled photon pairs + simple layout animation
% - Temperature tuning (affects n(λ,T) via linear thermo-optic approx)
% - Adjustable sample transmission spectrum (amplitude + gaussian notch)
% - Save MP4 button (writes animation to outputs/)
%
% IMPORTANT:
% - Replace the Sellmeier coefficients / poling period / L in params (or
%   create a params.mat in ./params/) to use exact numbers from your PDF.
% - Default coefficients are placeholders for demonstration only.

function PPLN
    % Create output folder if missing
    if ~exist('outputs','dir'); mkdir('outputs'); end
    if ~exist('params','dir'); mkdir('params'); end

    % ---- Default parameters (edit or save to params/params_custom.mat) ----
    P.lambda_p = 532e-9;        % pump center wavelength (m)
    P.lambda_s_target = 810e-9; % desired signal (m)
    P.pump_FWHM_nm = 1.0;       % pump FWHM (nm)
    P.L = 20e-3;                % crystal length (m)
    P.Lambda_poling = 7.0e-6;   % poling period (m) - Adjusted for sample generation
    P.Npairs = 2000;            % # Monte Carlo pairs
    P.Nvis = 400;               % # visible markers in animation
    P.w0 = 50e-6;               % pump waist (m)
    P.Aeff = 50e-12;            % effective modal area (m^2)
    P.d_eff = 2e-11;            % effective nonlinear (m/V)
    P.detect_s_eta = 0.6;       % signal detector efficiency
    P.detect_i_eta = 0.4;       % idler detector efficiency
    P.animate = true;
    P.saveMP4_default = true;   % default Save MP4 checkbox
    % Temperature (degC) - reference 25 C
    P.T0 = 25;
    P.T = 100;                  % Adjusted for sample generation

    % --- Sellmeier / n(λ,T) setup ---
    % You should replace the following placeholder coefficients with the
    % Sellmeier data from your PDF. I include a simple 3-term style as example.
    % Format: n^2(λ) = A + B/(λ^2 - C) + D/(λ^2 - E) + F/(λ^2 - G)
    Sellmeier.A = 5.35583; Sellmeier.B = 0.100473; Sellmeier.C = 0.20692^2;
    Sellmeier.D = 0.20692; Sellmeier.E = 0.0; Sellmeier.F = 0.0; Sellmeier.G = 0.0;
    Sellmeier.type = 'placeholder'; % set 'placeholder' until you replace with real set
    % Thermo-optic coefficient (approx dn/dT), one value per polarization:
    dn_dT = 3e-5; % 1/C (placeholder) - LiNbO3 typical ~ (1e-5 .. few e-5) range

    % ---- Geometry for animation ----
    geometry.zPPLN = 0.3;
    geometry.zDichroic = 0.36;
    geometry.zSample = 0.52;
    geometry.zDetSignal = 0.95;
    geometry.zDetIdler = 0.75;

    % UI: figure and axes
    fig = uifigure('Name','PPLN Designer v2','Position',[40 40 1300 760]);

    % Left control panel
    ctrlW = 330;
    pnl = uipanel(fig,'Title','Controls','Position',[10 10 ctrlW 760]);

    % Add controls (sliders + numeric)
    y0 = 720;
    uilabel(pnl,'Text','Poling period Λ (µm)','Position',[10 y0 180 18]);
    sLambda = uislider(pnl,'Limits',[1 30],'Value',P.Lambda_poling*1e6,'Position',[10 y0-30 300 3]);
    y0 = y0 - 65;
    uilabel(pnl,'Text','Crystal length L (mm)','Position',[10 y0 180 18]);
    sL = uislider(pnl,'Limits',[1 50],'Value',P.L*1e3,'Position',[10 y0-30 300 3]);
    y0 = y0 - 65;
    uilabel(pnl,'Text','Pump λ (nm)','Position',[10 y0 180 18]);
    sLambdaP = uislider(pnl,'Limits',[350 600],'Value',P.lambda_p*1e9,'Position',[10 y0-30 300 3]);
    y0 = y0 - 65;
    uilabel(pnl,'Text','Pump FWHM (nm)','Position',[10 y0 180 18]);
    sPW = uislider(pnl,'Limits',[0.1 10],'Value',P.pump_FWHM_nm,'Position',[10 y0-30 300 3]);
    y0 = y0 - 65;
    uilabel(pnl,'Text','Desired signal λ (nm)','Position',[10 y0 180 18]);
    sLS = uislider(pnl,'Limits',[400 1000],'Value',P.lambda_s_target*1e9,'Position',[10 y0-30 300 3]);
    y0 = y0 - 65;
    uilabel(pnl,'Text','# pairs (MC)','Position',[10 y0 180 18]);
    sNP = uislider(pnl,'Limits',[200 20000],'Value',P.Npairs,'Position',[10 y0-30 300 3]);
    sNP.MajorTicks = [200 1000 5000 10000 20000];
    y0 = y0 - 65;
    uilabel(pnl,'Text','Temperature (°C)','Position',[10 y0 180 18]);
    sT = uislider(pnl,'Limits',[-50 200],'Value',P.T,'Position',[10 y0-30 300 3]);
    y0 = y0 - 75;
    % Sample spectrum controls
    uilabel(pnl,'Text','Sample amplitude transmission (center)','Position',[10 y0 220 18]);
    edTs = uieditfield(pnl,'numeric','Value',0.7,'Position',[230 y0 70 22]);
    y0 = y0 - 45;
    uilabel(pnl,'Text','Sample phase shift (rad)','Position',[10 y0 220 18]);
    edPhi = uieditfield(pnl,'numeric','Value',0.0,'Position',[230 y0 70 22]);
    y0 = y0 - 45;
    uilabel(pnl,'Text','Sample spectral width (nm)','Position',[10 y0 220 18]);
    edSigmaSpec = uieditfield(pnl,'numeric','Value',50,'Position',[230 y0 70 22]);
    y0 = y0 - 55;
    % Save MP4 toggle and button
    cbSave = uicheckbox(pnl,'Text','Save MP4 on Run','Position',[10 y0 200 22],'Value',P.saveMP4_default);
    btnRun = uibutton(pnl,'push','Text','Resample & Animate','Position',[10 y0-55 140 34],'ButtonPushedFcn',@(btn,event) onRun());
    btnSaveSettings = uibutton(pnl,'push','Text','Save params to file','Position',[160 y0-55 150 34],'ButtonPushedFcn',@(btn,event) saveParams());
    y0 = y0 - 110;
    % Import Sellmeier button
    uilabel(pnl,'Text','Sellmeier / poling values:','Position',[10 y0 220 18]);
    btnLoadParams = uibutton(pnl,'push','Text','Load params/coeffs','Position',[10 y0-35 140 28],'ButtonPushedFcn',@(btn,event) loadParams());
    lblParams = uilabel(pnl,'Text','Using default placeholder Sellmeier. Click Load to override.','Position',[10 y0-70 300 36],'FontSize',10,'HorizontalAlignment','left');

    % Right/plots area
    axJSI = uiaxes(fig,'Position',[360 400 900 330]); title(axJSI,'JSI (approx)'); xlabel(axJSI,'\lambda_s (nm)'); ylabel(axJSI,'\lambda_i (nm)');
    axSlice = uiaxes(fig,'Position',[360 260 420 120]); title(axSlice,'Pump & phase-match (norm)'); xlabel(axSlice,'\lambda_s (nm)');
    axAnim = uiaxes(fig,'Position',[800 60 460 300]); axis(axAnim,'equal'); xlim(axAnim,[0 1]); ylim(axAnim,[-0.25 0.25]);
    title(axAnim,'Layout animation');

    % App state struct
    state.P = P; state.Sell = Sellmeier; state.dn_dT = dn_dT; state.geometry = geometry;
    state.controls = struct('sLambda',sLambda,'sL',sL,'sLambdaP',sLambdaP,'sPW',sPW,'sLS',sLS,'sNP',sNP,'sT',sT,'edTs',edTs,'edPhi',edPhi,'edSigmaSpec',edSigmaSpec);
    state.axes = struct('axJSI',axJSI,'axSlice',axSlice,'axAnim',axAnim);
    state.lblParams = lblParams; state.cbSave = cbSave; state.fig = fig;

    % connect slider callbacks for live update (but not too frequent)
    sLambda.ValueChangedFcn = @(s,e) quickUpdate();
    sL.ValueChangedFcn = @(s,e) quickUpdate();
    sLambdaP.ValueChangedFcn = @(s,e) quickUpdate();
    sPW.ValueChangedFcn = @(s,e) quickUpdate();
    sLS.ValueChangedFcn = @(s,e) quickUpdate();
    sNP.ValueChangedFcn = @(s,e) quickUpdate();
    sT.ValueChangedFcn = @(s,e) quickUpdate();

    % initial draw
    drawAll();

    % -------- nested functions --------
    function quickUpdate()
        % do light-weight redraw (JSI + slice) but not full resample
        readControlsToState();
        drawJSI_and_slice();
    end

    function onRun()
        readControlsToState();
        % draw and create animation (and possibly save mp4)
        drawAll();
        if cbSave.Value
            % save mp4 of animation
            fname = fullfile('outputs', sprintf('PPLN_anim_%s.mp4', datestr(now,'yyyymmdd_HHMMSS')));
            fprintf('Saving animation to %s (this may take a while)...\n', fname);
            writeAnimationToMP4(fname);
            uialert(fig, sprintf('Saved MP4 to:\n%s', fname), 'Saved');
        end
    end

    function saveParams()
        readControlsToState();
        paramsFile = fullfile('params','params_custom.mat');
        Psave = state.P; Sell = state.Sell; dn = state.dn_dT;
        save(paramsFile,'Psave','Sell','dn');
        uialert(fig,sprintf('Saved params to %s', paramsFile),'Saved');
    end

    function loadParams()
        [f,p] = uigetfile({'*.mat;*.txt;*.json','Parameter file (*.mat, *.txt, *.json)';'*.*','All files'},'Load params or sellmeier');
        if isequal(f,0); return; end
        pathf = fullfile(p,f);
        try
            ext = f(end-2:end);
            if strcmpi(ext,'mat')
                tmp = load(pathf);
                if isfield(tmp,'Psave'); state.P = tmp.Psave; end
                if isfield(tmp,'Sell'); state.Sell = tmp.Sell; end
                if isfield(tmp,'Sellmeier'); state.Sell = tmp.Sellmeier; end
                if isfield(tmp,'dn'); state.dn_dT = tmp.dn; end
                lblParams.Text = sprintf('Loaded params from %s', f);
            else
                % For text/json we would parse - for now show message
                lblParams.Text = 'Loaded file, but automatic parsing not implemented for txt/json. Use .mat with fields Psave,Sell.';
            end
            drawAll();
        catch ME
            uialert(fig,sprintf('Error loading params: %s', ME.message),'Load error');
        end
    end

    function readControlsToState()
        s = state.controls;
        state.P.Lambda_poling = s.sLambda.Value*1e-6;
        state.P.L = s.sL.Value*1e-3;
        state.P.lambda_p = s.sLambdaP.Value*1e-9;
        state.P.pump_FWHM_nm = s.sPW.Value;
        state.P.lambda_s_target = s.sLS.Value*1e-9;
        state.P.Npairs = round(s.sNP.Value);
        state.P.T = s.sT.Value;
        state.P.sample_Tcenter = s.edTs.Value;
        state.P.sample_phi = s.edPhi.Value;
        state.P.sample_sigma_nm = s.edSigmaSpec.Value;
    end

    function drawAll()
        readControlsToState();
        drawJSI_and_slice();
        drawAnimation();
    end

    function drawJSI_and_slice()
        P = state.P; Sell = state.Sell;
        ax = state.axes.axJSI; aslice = state.axes.axSlice;
        % build signal grid
        dw = 6e-9;
        lam_s_grid = linspace(P.lambda_s_target - dw, P.lambda_s_target + dw, 240);
        lam_i_grid = 1./(1./P.lambda_p - 1./lam_s_grid);
        % pump env
        sigma_lambda = (P.pump_FWHM_nm/2.355)*1e-9;
        pump_env = @(lam) exp(-(lam - P.lambda_p).^2/(2*sigma_lambda^2));
        % phase matching
        pm_amp = zeros(size(lam_s_grid));
        for ii = 1:length(lam_s_grid)
            ls = lam_s_grid(ii); li = lam_i_grid(ii);
            kp = wavevector(P.lambda_p, Sell, P.T, state.dn_dT);
            ks = wavevector(ls, Sell, P.T, state.dn_dT);
            ki = wavevector(li, Sell, P.T, state.dn_dT);
            Kq = 2*pi / P.Lambda_poling;
            Delta_k = kp - ks - ki - Kq;
            pm_amp(ii) = sinc_custom(Delta_k * P.L / 2);
        end
        JSI_curve = (abs(pm_amp).^2) .* (pump_env((lam_s_grid+lam_i_grid)/2)).^2;
        % show approximate 2D JSI (smeared)
        axes(ax); cla(ax);
        [LS, LI] = meshgrid(linspace(lam_s_grid(1),lam_s_grid(end),240), linspace(min(lam_i_grid),max(lam_i_grid),240));
        % interpolate JSI_curve to LS
        Jvec = interp1(lam_s_grid, JSI_curve, LS(1,:),'linear',0);
        imagesc(ax, (LS(1,:)*1e9), (LI(:,1)*1e9), repmat(Jvec, size(LI,1),1));
        axis(ax,'xy'); xlabel(ax,'\lambda_s (nm)'); ylabel(ax,'\lambda_i (nm)');
        title(ax,'Approx JSI (along energy-conservation curve)');
        hold(ax,'on'); plot(ax, P.lambda_s_target*1e9, (1/(1/P.lambda_p - 1/P.lambda_s_target))*1e9, 'r+','MarkerSize',10); hold(ax,'off');

        % slice: normalized pump and PM
        axes(aslice); cla(aslice);
        plot(aslice, lam_s_grid*1e9, pump_env((lam_s_grid+lam_i_grid)/2)/max(pump_env((lam_s_grid+lam_i_grid)/2)),'--b','DisplayName','Pump env (norm)'); hold(aslice,'on');
        plot(aslice, lam_s_grid*1e9, (abs(pm_amp).^2)/max(abs(pm_amp).^2), '-r','DisplayName','PM (norm)');
        xlabel(aslice,'\lambda_s (nm)'); ylabel(aslice,'Normalized'); legend(aslice,'Location','best'); grid(aslice,'on');
    end

    function drawAnimation()
        P = state.P; Sell = state.Sell;
        % build sampling grid
        dw = 6e-9;
        lam_s_grid = linspace(P.lambda_s_target - dw, P.lambda_s_target + dw, 600);
        lam_i_grid = 1./(1./P.lambda_p - 1./lam_s_grid);
        sigma_lambda = (P.pump_FWHM_nm/2.355)*1e-9;
        pump_env = @(lam) exp(-(lam - P.lambda_p).^2/(2*sigma_lambda^2));
        pm_amp = zeros(size(lam_s_grid));
        for ii = 1:length(lam_s_grid)
            ls = lam_s_grid(ii); li = lam_i_grid(ii);
            kp = wavevector(P.lambda_p, Sell, P.T, state.dn_dT);
            ks = wavevector(ls, Sell, P.T, state.dn_dT);
            ki = wavevector(li, Sell, P.T, state.dn_dT);
            Kq = 2*pi / P.Lambda_poling;
            Delta_k = kp - ks - ki - Kq;
            pm_amp(ii) = sinc_custom(Delta_k * P.L / 2);
        end
        JSI_curve = (abs(pm_amp).^2) .* (pump_env((lam_s_grid+lam_i_grid)/2)).^2;
        JSI_curve(JSI_curve<0)=0;
        if sum(JSI_curve) == 0
            uialert(state.fig, 'No photon pairs generated with current parameters. Adjust settings.', 'Warning');
            return;
        end
        pdf = JSI_curve / trapz(lam_s_grid, JSI_curve);
        cdf = cumtrapz(lam_s_grid, pdf); cdf = cdf / cdf(end);

        % sample
        N = P.Npairs;
        r = rand(1,N);
        lam_s_samples = interp1(cdf, lam_s_grid, r, 'linear', 'extrap');
        lam_i_samples = 1./(1./P.lambda_p - 1./lam_s_samples);

        % angular emission approx
        sigma_theta = 0.5*(P.lambda_p/(pi*P.w0));
        theta_s = sigma_theta .* randn(1,N);
        theta_i = -theta_s + 1e-3*randn(1,N);

        % sample structure
        samples.ls = lam_s_samples; samples.li = lam_i_samples;
        samples.theta_s = theta_s; samples.theta_i = theta_i;
        samples.alive_i = true(1,N); samples.alive_s = true(1,N);
        samples.applied = false(1,N);

        % animation settings
        ax = state.axes.axAnim;
        cla(ax); hold(ax,'on');
        % static layout
        plot(ax, [geometry.zPPLN geometry.zPPLN], [-0.18 0.18], 'k','LineWidth',3); text(ax, geometry.zPPLN, 0.2, 'PPLN','HorizontalAlignment','center');
        plot(ax, [geometry.zDichroic geometry.zDichroic], [-0.18 0.18], 'm--','LineWidth',2); text(ax, geometry.zDichroic, 0.2,'Dichroic','HorizontalAlignment','center');
        plot(ax, [geometry.zSample geometry.zSample], [-0.18 0.18], 'r','LineWidth',2); text(ax, geometry.zSample, 0.2,'Sample','HorizontalAlignment','center');
        plot(ax, geometry.zDetSignal, 0.0, 'gs','MarkerFaceColor','g'); text(ax, geometry.zDetSignal,0.05,'Signal Det','HorizontalAlignment','center');
        plot(ax, geometry.zDetIdler, -0.12, 'ys','MarkerFaceColor','y'); text(ax, geometry.zDetIdler,-0.15,'Idler Det','HorizontalAlignment','center');

        % animate subset of pairs for clarity
        idxVis = randperm(N, min(P.Nvis, N));
        nframes = 200;
        for f = 1:nframes
            cla(ax); hold(ax,'on');
            % redraw static elements
            plot(ax, [geometry.zPPLN geometry.zPPLN], [-0.18 0.18], 'k','LineWidth',3);
            plot(ax, [geometry.zDichroic geometry.zDichroic], [-0.18 0.18], 'm--','LineWidth',2);
            plot(ax, [geometry.zSample geometry.zSample], [-0.18 0.18], 'r','LineWidth',2);
            plot(ax, geometry.zDetSignal, 0.0, 'gs','MarkerFaceColor','g');
            plot(ax, geometry.zDetIdler, -0.12, 'ys','MarkerFaceColor','y');

            tnorm = (f-1)/(nframes-1);
            for ii = 1:length(idxVis)
                kidx = idxVis(ii);
                if tnorm <= 0.35
                    % between PPLN and dichroic
                    xs = geometry.zPPLN + (geometry.zDichroic - geometry.zPPLN) * (tnorm/0.35);
                    ys = +0.01*(tnorm/0.35);
                    xi = xs; yi = -0.01*(tnorm/0.35);
                else
                    % arms after dichroic
                    frac = min((tnorm-0.35)/0.65,1);
                    xs = geometry.zDichroic + (geometry.zDetSignal - geometry.zDichroic)*frac;
                    ys = 0.05*(1-frac);
                    xi = geometry.zDichroic + (geometry.zDetIdler - geometry.zDichroic)*frac;
                    yi = -0.02 - 0.1*frac;
                    % when idler crosses sample_x, apply sample transmission
                    if ~isfield(samples,'applied') || ~samples.applied(kidx)
                        % compute when xi passes sample location (approx)
                        fracSample = (geometry.zSample - geometry.zDichroic) / (geometry.zDetIdler - geometry.zDichroic);
                        if frac >= fracSample
                            % wavelength-dependent sample transmission
                            Ts = sampleTransmission(samples.li(kidx), state.P);
                            if rand > Ts
                                samples.alive_i(kidx) = false;
                            end
                            samples.applied(kidx) = true;
                        end
                    end
                end
                if samples.alive_s(kidx)
                    plot(ax, xs, ys, 'bo','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');
                end
                if samples.alive_i(kidx)
                    plot(ax, xi, yi, 'ro','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');
                else
                    % show lost idler near sample
                    plot(ax, geometry.zSample, -0.06, 'x','Color',[0.6 0.6 0.6]);
                end
            end
            xlim(ax,[0 1]); ylim(ax,[-0.25 0.25]); xlabel(ax,'x (m)'); ylabel(ax,'y (m)');
            drawnow;
        end
    end

    function writeAnimationToMP4(fname)
        % Create a video writer of the same animation above but capture frames
        P = state.P; Sell = state.Sell;
        v = VideoWriter(fname,'MPEG-4'); v.FrameRate = 30; open(v);
        % call drawAnimation but capture frames
        % Re-run Monte Carlo sampling
        dw = 6e-9;
        lam_s_grid = linspace(P.lambda_s_target - dw, P.lambda_s_target + dw, 600);
        lam_i_grid = 1./(1./P.lambda_p - 1./lam_s_grid);
        sigma_lambda = (P.pump_FWHM_nm/2.355)*1e-9;
        pump_env = @(lam) exp(-(lam - P.lambda_p).^2/(2*sigma_lambda^2));
        pm_amp = zeros(size(lam_s_grid));
        for ii = 1:length(lam_s_grid)
            ls = lam_s_grid(ii); li = lam_i_grid(ii);
            kp = wavevector(P.lambda_p, Sell, P.T, state.dn_dT);
            ks = wavevector(ls, Sell, P.T, state.dn_dT);
            ki = wavevector(li, Sell, P.T, state.dn_dT);
            Kq = 2*pi / P.Lambda_poling;
            Delta_k = kp - ks - ki - Kq;
            pm_amp(ii) = sinc_custom(Delta_k * P.L / 2);
        end
        JSI_curve = (abs(pm_amp).^2) .* (pump_env((lam_s_grid+lam_i_grid)/2)).^2;
        pdf = JSI_curve; pdf(pdf<0)=0;
        if sum(pdf) == 0
            uialert(state.fig, 'No photon pairs generated with current parameters. Adjust settings.', 'Warning');
            close(v);
            return;
        end
        pdf = pdf / trapz(lam_s_grid, pdf);
        cdf = cumtrapz(lam_s_grid, pdf); cdf = cdf / cdf(end);
        N = P.Npairs;
        r = rand(1,N);
        lam_s_samples = interp1(cdf, lam_s_grid, r, 'linear', 'extrap');
        lam_i_samples = 1./(1./P.lambda_p - 1./lam_s_samples);
        sigma_theta = 0.5*(P.lambda_p/(pi*P.w0));
        theta_s = sigma_theta .* randn(1,N);
        theta_i = -theta_s + 1e-3*randn(1,N);
        samples.ls = lam_s_samples; samples.li = lam_i_samples;
        samples.theta_s = theta_s; samples.theta_i = theta_i;
        samples.alive_i = true(1,N); samples.alive_s = true(1,N);
        samples.applied = false(1,N);
        idxVis = randperm(N, min(P.Nvis, N));
        % frames
        nframes = 300;
        for f = 1:nframes
            % create the same plot on an invisible figure and capture
            fh = figure('Visible','off','Position',[100 100 900 600]);
            ax = axes(fh); hold(ax,'on');
            % static
            plot(ax, [geometry.zPPLN geometry.zPPLN], [-0.18 0.18], 'k','LineWidth',3);
            plot(ax, [geometry.zDichroic geometry.zDichroic], [-0.18 0.18], 'm--','LineWidth',2);
            plot(ax, [geometry.zSample geometry.zSample], [-0.18 0.18], 'r','LineWidth',2);
            plot(ax, geometry.zDetSignal, 0.0, 'gs','MarkerFaceColor','g');
            plot(ax, geometry.zDetIdler, -0.12, 'ys','MarkerFaceColor','y');
            % draw subset
            tnorm = (f-1)/(nframes-1);
            for ii = 1:length(idxVis)
                kidx = idxVis(ii);
                if tnorm <= 0.35
                    xs = geometry.zPPLN + (geometry.zDichroic - geometry.zPPLN) * (tnorm/0.35);
                    ys = +0.01*(tnorm/0.35);
                    xi = xs; yi = -0.01*(tnorm/0.35);
                else
                    frac = min((tnorm-0.35)/0.65,1);
                    xs = geometry.zDichroic + (geometry.zDetSignal - geometry.zDichroic)*frac;
                    ys = 0.05*(1-frac);
                    xi = geometry.zDichroic + (geometry.zDetIdler - geometry.zDichroic)*frac;
                    yi = -0.02 - 0.1*frac;
                    if ~isfield(samples,'applied') || ~samples.applied(kidx)
                        fracSample = (geometry.zSample - geometry.zDichroic) / (geometry.zDetIdler - geometry.zDichroic);
                        if frac >= fracSample
                            Ts = sampleTransmission(samples.li(kidx), state.P);
                            if rand > Ts
                                samples.alive_i(kidx) = false;
                            end
                            samples.applied(kidx) = true;
                        end
                    end
                end
                if samples.alive_s(kidx)
                    plot(ax, xs, ys, 'bo','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');
                end
                if samples.alive_i(kidx)
                    plot(ax, xi, yi, 'ro','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');
                else
                    plot(ax, geometry.zSample, -0.06, 'x','Color',[0.6 0.6 0.6]);
                end
            end
            xlim(ax,[0 1]); ylim(ax,[-0.25 0.25]); xlabel(ax,'x (m)'); ylabel(ax,'y (m)');
            frame = getframe(fh);
            writeVideo(v, frame);
            close(fh);
        end
        close(v);
    end

    function k = wavevector(lambda, Sell, T, dn_dT)
        % Compute wavevector k = 2*pi * n(lambda,T) / lambda
        nlam = sellmeier_n(lambda, Sell);
        % thermo-optic linear approx
        nlam = nlam + dn_dT * (T - state.P.T0);
        k = 2*pi * nlam ./ lambda;
    end

    function Ts = sampleTransmission(lambda_i, Pparams)
        % wavelength-dependent sample transmission modeled as Gaussian around
        % center with amplitude set by sample_Tcenter and width sample_sigma_nm
        lam = lambda_i*1e9;
        center = Pparams.lambda_s_target*1e9*1.915; % example center roughly; user can adjust
        sigma = Pparams.sample_sigma_nm;
        % Gaussian line shape (amplitude)
        Ts0 = Pparams.sample_Tcenter;
        Ts = Ts0 + (1-Ts0)*exp(-0.5*((lam-center)/sigma).^2);
        Ts = min(max(Ts,0),1);
    end

    function y = sinc_custom(x)
        y = ones(size(x));
        nz = x~=0;
        y(nz) = sin(x(nz))./x(nz);
    end

    function nv = sellmeier_n(lambda, Sell)
        % placeholder Sellmeier: n^2 = A + B/(lambda^2 - C) + ...
        % expect lambda in meters. Sell fields A,B,C,D,E,F,G optional.
        lam2 = (lambda*1e6).^2; % convert m -> microns^2 for typical Sellmeier units if coefficients fitted that way
        % NOTE: check units of your Sellmeier set. Adjust conversions accordingly.
        % We'll assume the placeholder coefficients were given for microns.
        nv2 = Sell.A;
        if isfield(Sell,'B') && isfield(Sell,'C')
            nv2 = nv2 + Sell.B./(lam2 - Sell.C);
        end
        if isfield(Sell,'D') && isfield(Sell,'E') && Sell.D~=0
            nv2 = nv2 + Sell.D./(lam2 - Sell.E);
        end
        nv = sqrt(abs(nv2));
    end
end
