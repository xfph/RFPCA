function [bp, X, opts] = mvt_ini(X, opts)
% initialization of parameters

tic;
[d(1), d(2), N]=size(X);
bp.d = d; 
bp.nu = 10;

if (~exist('opts', 'var')) opts = []; end
if ~isfield(opts, 'center') opts.center = 0; end
if ~isfield(opts, 'alg') opts.alg = 'PX-ECME'; end
if ~isfield(opts, 'ini') opts.ini = 'random'; end
if ~isfield(opts, 'tol') opts.tol = 1e-10; end
if ~isfield(opts, 'regv') opts.regv = 1e-6; end

bp.M = mean(X, 3);
if opts.center
    X = X - bp.M; 
end

if strcmp(opts.ini,'2dsvd') || strcmp(opts.ini,'adhoc')
    X1 = X(:, :); F{1} = X1 * X1';
    X2 = permute(X, [2, 1, 3]); X2 = X2(:, :); F{2} = X2 * X2';
end

%rng(1,'twister');
switch opts.ini
    case 'random' 
        Sc = rand(d(1)); bp.Sc = Sc*Sc';
        Sr = rand(d(2)); bp.Sr = Sr*Sr';
    case '2dsvd'
        bp.Sc=F{1}/(N*dr);
        bp.Sr=F{2}/(N*dc);
    case 'adhoc'
        bp.Sc=F{1};
        bp.Sr=F{2};
end
opts.time.ini = toc;       
        