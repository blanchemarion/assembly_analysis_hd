function jitter=assembly_deterministic_jitter(n,width)
%ASSEMBLY_DETERMINISTIC_JITTER Symmetric reproducible offsets for dot plots.
if n<=1,jitter=zeros(n,1);else,jitter=linspace(-width,width,n)';end
end
