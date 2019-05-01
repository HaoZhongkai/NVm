function S = Cal_Square(data,Px,W)
%%%data为数据，Px为在当前参数下拟合出的值
S = (data-Px)'*diag(W)*(data-Px);