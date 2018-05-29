function [w] = FKL_weights(Kernels_list,adjmat,dim,r_lamda)
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
% [1]He J, Chang S F, Xie L. Fast kernel learning for spatial pyramid matching[C]// 
%             Computer Vision and Pattern Recognition, 2008. CVPR 2008. IEEE Conference on. IEEE, 2008:1-7.
% adjmat : binary adjacency matrix
% dim    : dimension (1 - rows, 2 - cols)
% r_lamda: Regularized item (20000)
num_kernels = size(Kernels_list,3);
%r_lamda = 20000;
weight_v = zeros(num_kernels,1);

y = adjmat;
    % Graph based kernel
if dim == 1
        ga = y*y';
else
        ga = y'*y;
end

%Knormalized Kernels

%for i=1:num_kernels
	
%	S=Knormalized(Kernels_list(:,:,i));
%	Kernels_list(:,:,i) = S;

%end
%N_U = size(y,dim);
%l=ones(N_U,1);
%U = eye(N_U) - (l*l')/N_U;

M = zeros(num_kernels,num_kernels);

for i=1:num_kernels
	for j=1:num_kernels
		kk1 = Kernels_list(:,:,i);
		kk2 = Kernels_list(:,:,j);
		mm = trace(kk1'*kk2);
		M(i,j) = mm;
	end
end

a = zeros(num_kernels,1);

for i=1:num_kernels

	kk = Kernels_list(:,:,i);
	aa = trace(ga'*kk);
	a(i) = aa;
end

v = randn(num_kernels,1);
M_A = (M + r_lamda*eye(num_kernels));
cvx_begin
    variable v(num_kernels,1);
    minimize( v'*M_A*v - 2*a'*v );
	v >= 0;
	sum(v)==1;
cvx_end

%w = v /sum(v);
%w = v /norm(v,2);
w = v;