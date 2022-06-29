function [bp, x, opts] = mt_ini(x, opts)
% initialization of parameters

tic;
[d, N]=size(x);
bp.d = d; 
bp.nu = 4;
bp.tau = rand(1,N);

if (~exist('opts', 'var')) opts = []; end
if ~isfield(opts, 'center') opts.center = 0; end
if ~isfield(opts, 'alg') opts.alg = 'ECME'; end
if ~isfield(opts, 'ini') opts.ini = 'random'; end
if ~isfield(opts, 'tol') opts.tol = 1e-10; end
if ~isfield(opts, 'regv') opts.regv = 1e-6; end

bp.M = mean(x, 2);
if opts.center
    x = x - bp.M; 
end

if  strcmp(opts.ini,'adhoc')
    F = x * x';
end

switch opts.ini
    case 'random' 
        S = rand(d); bp.S = S*S';
    case 'adhoc'
        bp.S = F;
end
opts.time.ini = toc;       
        