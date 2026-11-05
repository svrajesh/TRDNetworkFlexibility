function C = regularizeCov(C)

% Symmetrize

C = (C+C')/2;


% Scale-dependent ridge regularization

reg = 1e-6 * trace(C)/size(C,1);


if ~isfinite(reg) || reg==0
    reg = 1e-8;
end


C = C + reg*eye(size(C));

end
