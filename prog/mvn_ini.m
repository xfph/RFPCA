function [bp, X, opts] = mvn_ini(X, opts)
% initialization of parameters

tic;
[d(1), d(2), N]=size(X);
bp.d = d; 

if (~exist('opts', 'var')) opts = []; end
if ~isfield(opts, 'center') opts.center = 1; end
if ~isfield(opts, 'alg') opts.alg = 'CM'; end
if ~isfield(opts, 'ini') opts.ini = 'random'; end
if ~isfield(opts, 'tol') opts.tol = 1e-10; end
if ~isfield(opts, 'regv') opts.regv = 1e-6; end

bp.M = mean(X, 3);
if opts.center
    X = X - bp.M; 
end

if strcmp(opts.ini,'2dsvd') || strcmp(opts.ini,'adhoc')
    X2 = permute(X, [2, 1, 3]); X2 = X2(:,:);
    F{2} = X2 * X2';
end

%rng(1,'twister');
switch opts.ini
    case 'random' 
        Sr = rand(d(2)); bp.Sr = Sr*Sr';
    case '2dsvd'
        bp.Sr=F{2}/(N*dc);
    case 'adhoc'
        bp.Sr=F{2};
end
opts.time.ini = toc;       
        