function S = Cal_Square(data,Px,W)
%%%dataΪ���ݣ�PxΪ�ڵ�ǰ��������ϳ���ֵ
S = (data-Px)'*diag(W)*(data-Px);