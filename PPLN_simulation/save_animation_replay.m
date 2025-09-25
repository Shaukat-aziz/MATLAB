% save_animation_replay.m
% Robust MP4 writer for PPLN animation.
% Usage: run save_animation_replay from project folder after running simulation_driver
% It will load outputs/R_vs_dli.mat if present; otherwise you can pass P,dli_vec,Rvec manually.

function save_animation_replay()
    try
        if ~exist('outputs','dir'), mkdir('outputs'); end
        datafile = fullfile('outputs','R_vs_dli.mat');
        if ~exist(datafile,'file')
            error('Missing data file: %s. Run simulation_driver_quick (no animation) first.', datafile);
        end
        S = load(datafile); dli_vec = S.dli_vec; Rvec = S.Rvec; P = S.P;
    catch ME
        error('Could not load saved data: %s\n', ME.message);
    end

    % prepare video file
    outname = fullfile('outputs', sprintf('PPLN_simulation_fixed_%s.mp4', datestr(now,'yyyymmdd_HHMMSS')));
    v = VideoWriter(outname, 'MPEG-4');
    v.FrameRate = 25;
    open(v);

    % We'll render frames into an offscreen (invisible) figure to avoid desktop issues
    fh = figure('Visible','off','Position',[100 100 1200 600],'Renderer','painters'); % try software renderer
    ax1 = subplot(1,2,1,'Parent',fh); axis(ax1,'equal'); xlim(ax1,[0 1]); ylim(ax1,[-0.3 0.3]); hold(ax1,'on');
    ax2 = subplot(1,2,2,'Parent',fh); hold(ax2,'on'); xlabel(ax2,'\Delta l_i (mm)'); ylabel(ax2,'R (norm)');
    plot(ax2, dli_vec*1e3, Rvec, 'Color',[0.6 0.6 0.6]);
    hh = plot(ax2, dli_vec(1)*1e3, Rvec(1), 'ro','MarkerFaceColor','r');

    % static layout coordinates (same as animate_beams)
    zPPLN = 0.3; zD = 0.36; zS = 0.52; zDetS = 0.95; zDetI = 0.75;
    % Precompute number of frames
    nframes = max(100, P.animationFrames);
    Npairs_vis = 250; % how many points to draw each frame (keeps file smallish)

    fprintf('Writing video to %s with %d frames... \n', outname, nframes);
    for fr = 1:nframes
        % draw left axes
        cla(ax1); hold(ax1,'on');
        plot(ax1, [zPPLN zPPLN], [-0.25 0.25], 'k','LineWidth',3);
        plot(ax1, [zD zD], [-0.25 0.25], 'm--','LineWidth',2);
        plot(ax1, [zS zS], [-0.25 0.25], 'r','LineWidth',2);
        plot(ax1, zDetS, 0.0, 'gs','MarkerFaceColor','g');
        plot(ax1, zDetI, -0.12, 'ys','MarkerFaceColor','y');

        % draw some moving particles (pseudo-random positions based on frame)
        rng(100 + fr); % deterministic per-frame seed for reproducibility
        tnorm = (fr-1)/(nframes-1);
        for k = 1:Npairs_vis
            r = rand();
            if tnorm <= 0.35
                xs = zPPLN + (zD - zPPLN) * (tnorm / 0.35) * (0.8 + 0.4*randn()*0.02);
                ys = 0.02 * (tnorm / 0.35) * (1 + 0.2*(rand-0.5));
                xi = xs; yi = -0.02 * (tnorm / 0.35) * (1 + 0.2*(rand-0.5));
            else
                frac = min((tnorm - 0.35)/0.65, 1);
                xs = zD + (zDetS - zD) * frac * (0.98 + 0.04*randn()*0.02);
                ys = 0.05 * (1 - frac) * (1 + 0.2*(rand-0.5));
                xi = zD + (zDetI - zD) * frac * (0.98 + 0.04*randn()*0.02);
                yi = -0.02 - 0.12 * frac + 0.01*(rand-0.5);
            end
            plot(ax1, xs, ys, 'bo', 'MarkerFaceColor','b', 'MarkerSize',5);
            plot(ax1, xi, yi, 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
        end
        title(ax1, sprintf('Layout (frame %d/%d)  \\Delta l_i = %.3f mm', fr, nframes, dli_vec(max(1, round((fr-1)/(nframes-1)*(length(dli_vec)-1))+1))*1e3));
        xlim(ax1,[0 1]); ylim(ax1,[-0.3 0.3]);

        % update cursor on R plot
        cur_idx = max(1, round( (fr-1)/(nframes-1)*(length(dli_vec)-1) ) + 1);
        set(hh, 'XData', dli_vec(cur_idx)*1e3, 'YData', Rvec(cur_idx));

        % capture frame from invisible figure
        drawnow; % force render
        frame = getframe(fh);
        try
            writeVideo(v, frame);
        catch WVerr
            warning('writeVideo failed at frame %d: %s', fr, WVerr.message);
            break;
        end
    end

    close(v);
    close(fh);
    fprintf('Finished writing video: %s\n', outname);
end
