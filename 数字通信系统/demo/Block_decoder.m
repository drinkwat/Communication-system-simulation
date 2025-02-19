function est_code = Block_decoder(n,k,rec,H,d)
%输入
%n：码字长度
%k：信息组长度
%rec：译码器接收的磧判决码字
%H：分组码的校验矩阵
%d：分组码的最小距离
%输出
%est_code：译码输出码字

    [M,N] = size (H);
    if N~=n
       disp('Parameter of Code Length is error.\n');
    end
    if M~=n-k
       disp('Parameter of Parity Length is error.\n');
    end
    t=fix((d-1)/2);%纠正t个随机错误
    %初始化并计算不同个数H矩阵列向量的所有组合的个数
    num = zeros(1,t);
    for idx = 1:t
       num(idx) = factorial(n)/factorial(n-idx)/factorial(idx);
    end%计算纠正1~t个错误的需要H矩阵列向量的组合数

    maxnum = max(num);  %计算最大组合数
    out = zeros(t,maxnum,t);%分配空间;第1位代表纠几位，第2位代表错误图样个数，第3位代表每个错误图样里的错误
    out(1,1:num(1),1) = 1:n;%单个H列向量索引
    out = compound(out,n,num,t,1);%2到t个H列向量组合的所有可能索引集合
    
    est_code = rec;%如果超出纠错能力则不进行纠错，原样输出

    Hcom = zeros(1,n-k);    %初始化H组合列向量
    E = zeros(1,N); %初始化错误图样
    S = rem(rec*H',2);  %计算伴随式
    if find(S)  %伴随式不为零，表示码子中有错
        for err=1:t %在分组码的纠错能力范围内，查找与伴随式值匹配的H列向量或器组合。
            %从1位错误开始纠错
           for idx=1:num(err) %逐个检查H列向量组合。检查所有的错误图样
               Hcom = zeros(1, n-k);
               for j = 1:err %计算H组合列向量；针对具体的错误图样包含err个错码
                  Hcom = rem(Hcom + H(:,out(err,idx,j))',2);
               end
               %Hcom为对应H的列加和
               if(sum(rem(S+Hcom, 2)) == 0) %找到匹配的H列向量组合
                  E(out(err,idx,1:err)) = 1;    %估计错误图样，设置E对应的列为1
                  est_code = rem(rec+E, 2); %纠错
                  return;
               end
           end
        end
    else %无错，直接输出
        est_code = rec;
    end
end