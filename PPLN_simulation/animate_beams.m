% animate_beams.m
% Visual animation combining a 2D layout (pump -> PPLN -> split arms -> sample -> detectors)
% and overlay of R(dli) evolving along frames.
%
% Usage: animate_beams(P, dli_vec, Rvec)

function animate_beams(P, dli_vec, Rvec)
    % create figure
    fh = figure('Name','PPLN Animation','NumberTitle','off','Position',[100 100 1200 600]);
    ax1 = subplot(1,2,1); axis(ax1,'equal'); xlim(ax1,[0 1]); ylim(ax1,[-0.3 0.3]); hold(ax1,'on');
    title(ax1,'Optical Layout (cartoon)');
    ax2 = subplot(1,2,2); hold(ax2,'on'); xlabel(ax2,'\Delta l_i (mm)'); ylabel(ax2,'R (norm)');
    plot(ax2, dli_vec*1e3, Rvec, 'Color',[0.6 0.6 0.6]);
    hh = plot(ax2, dli_vec(1)*1e3, Rvec(1), 'ro','MarkerFaceColor','r');

    % static elements in layout
    zPPLN = 0.3; zD = 0.36; zS = 0.52; zDetS = 0.95; zDetI = 0.75;
    plot(ax1, [zPPLN zPPLN], [-0.25 0.25], 'k','LineWidth',3); text(ax1,zPPLN,0.27,'PPLN','HorizontalAlignment','center');
    plot(ax1, [zD zD], [-0.25 0.25], 'm--','LineWidth',2); text(ax1,zD,0.27,'Dichroic','HorizontalAlignment','center');
    plot(ax1, [zS zS], [-0.25 0.25], 'r','LineWidth',2); text(ax1,zS,0.27,'Sample','HorizontalAlignment','center');
    plot(ax1, zDetS, 0.0, 'gs','MarkerFaceColor','g'); text(ax1, zDetS, 0.05,'Signal Det','HorizontalAlignment','center');
    plot(ax1, zDetI, -0.12, 'ys','MarkerFaceColor','y'); text(ax1, zDetI, -0.16,'Idler Det','HorizontalAlignment','center');

    nframes = P.animationFrames;
    Npairs = 300; % visual pairs
    idx = randperm(5000, Npairs); % random seeds for visualization

    for fr = 1:nframes
        cla(ax1); hold(ax1,'on');
        % redraw static
        plot(ax1, [zPPLN zPPLN], [-0.25 0.25], 'k','LineWidth',3);
        plot(ax1, [zD zD], [-0.25 0.25], 'm--','LineWidth',2);
        plot(ax1, [zS zS], [-0.25 0.25], 'r','LineWidth',2);
        plot(ax1, zDetS, 0.0, 'gs','MarkerFaceColor','g');
        plot(ax1, zDetI, -0.12, 'ys','MarkerFaceColor','y');

        % animate a handful of rays moving from PPLN to detectors
        tnorm = (fr-1)/(nframes-1);
        for k = 1:round(Npairs/15)
            % parametric positions: before dichroic or after
            if tnorm <= 0.35
                xs = zPPLN + (zD - zPPLN) * (tnorm / 0.35);
                ys = 0.02 * (tnorm / 0.35);
                xi = xs; yi = -0.02 * (tnorm / 0.35);
            else
                frac = min((tnorm - 0.35)/0.65, 1);
                xs = zD + (zDetS - zD) * frac;
                ys = 0.05 * (1 - frac);
                xi = zD + (zDetI - zD) * frac;
                yi = -0.02 - 0.12 * frac;
            end
            plot(ax1, xs, ys, 'bo', 'MarkerFaceColor','b');
            plot(ax1, xi, yi, 'ro', 'MarkerFaceColor','r');
        end

        % update R overlay cursor
        cur_idx = max(1, round( (fr-1)/(nframes-1)*(length(dli_vec)-1) ) + 1);
        set(hh, 'XData', dli_vec(cur_idx)*1e3, 'YData', Rvec(cur_idx));

        title(ax1, sprintf('Layout (frame %d/%d)  \\Delta l_i = %.2f mm', fr, nframes, dli_vec(cur_idx)*1e3));
        drawnow;
    end

    % optionally save video
    try
        if ~isempty(P.videoFile)
            fprintf('Saving animation to %s ...\n', P.videoFile);
            v = VideoWriter(P.videoFile,'MPEG-4'); v.FrameRate = 25; open(v);
            for fr = 1:nframes
                % reproduce frames for video (cheap re-render)
                % to keep code simple we re-copy current figure to video
                frame = getframe(fh);
                writeVideo(v, frame);
            end
            close(v);
            fprintf('Saved video.\n');
        end
    catch ME
        warning('Could not save video: %s', ME.message);
    end
end
