% export_video.m
% Utility to re-render saved frames or re-run animation and produce a video.
% Simple wrapper that calls animate_beams with file saving enabled inside animate_beams.
function export_video()
    if ~exist('outputs','dir'); mkdir('outputs'); end
    P = params();
    load(fullfile('outputs','R_vs_dli.mat'),'dli_vec','Rvec');
    if isempty(P.videoFile)
        P.videoFile = fullfile('outputs', ['PPLN_anim_' datestr(now,'yyyymmdd_HHMMSS') '.mp4']);
    end
    animate_beams(P, dli_vec, Rvec);
    fprintf('export_video: done.\n');
end
