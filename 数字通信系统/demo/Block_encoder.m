function code = Block_encoder(n,k,msg,G)
%输入
%n：码字长度
%k：信息组长度
%msg,编码输入信息组
%输岀
%code：編码输出码字

    [M,N] = size(G);%获取生成矩阵的维数
    if N~=n
       disp('Parameter of Code Length is error.\n');
    end
    if M~=k
       disp('Parameter of Info. Length is error.\n');
    end
    %计算編码输出
    code = rem(msg*G,2);%模2
    
end