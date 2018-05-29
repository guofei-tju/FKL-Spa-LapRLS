clear;
load('You_dataset.mat');
MU = 0.005;GAMMA=0.002;
y = miRNA_disease_Y;

K1 = [];
K1(:,:,1)=miRNA_Function_S;
K1(:,:,2)=miRNA_Sequences_Needle_S;
K2 = [];
K2(:,:,1)=disease_Function_S;
K2(:,:,2)=disease_Sem_S;


y = miRNA_disease_Y;
nfolds =5; nruns=1;
crossval_idx = crossvalind('Kfold',y(:),nfolds);
for fold = 1:nfolds
	y_train = miRNA_disease_Y;
	test_idx  = find(crossval_idx==fold);
	y_train(test_idx) = 0;
	K1(:,:,3)=kernel_gip(y_train,1, 1);
	K2(:,:,3)=kernel_gip(y_train,2, 1);
 
	[weight_v1] = FKL_weights(K1,y_train,1,20000);
	K_COM1 = combine_kernels(weight_v1, K1);		

	[weight_v2] = FKL_weights(K2,y_train,2,20000);
	K_COM2 = combine_kernels(weight_v2, K2);
	
	[F_1] = LapRLS(K_COM1,K_COM2,y_train, 2^(-4),40,1);
	
	y(test_idx)= F_1(test_idx);
end

[X_1,Y_1,tpr,aupr_F_1] = perfcurve(miRNA_disease_Y(:),y(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_F_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(miRNA_disease_Y(:),y(:),1);
