function out = phantomBlank(data, fs, f, wide)
% Apply the same pre/post/context excision + pchip on a regular pulse grid
% with random phase, to OFF data (no pulses present).
pre = 25; post = 50 + 50*wide; ctx = 20;
T = size(data,2); out = data;
period = fs/f;
locs = round(randi(round(period)) : period : T);
locs = locs(locs > pre+ctx & locs < T-post-ctx);
for k = 1:numel(locs)
    gap   = (locs(k)-pre):(locs(k)+post);
    known = [(gap(1)-ctx):(gap(1)-1), (gap(end)+1):(gap(end)+ctx)];
    out(:,gap) = pchip(known, data(:,known), gap);
end
end