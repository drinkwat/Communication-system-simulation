%计算伴随式组合函数
function out = compound(out, n, num, d, idx)
% d,纠错最大数
% n,码长度
% num, 存储了所有可能的组合数
    if idx == d
        return;
    end

    flag= 1;
    cnt = 1;
    for i=1:num(idx) %组合数
        for j=out(idx,i,idx)+1:n
           for k=1:idx %纠错数目
               if(j==out(idx, i, idx))
                   flag =0;
               end
           end
           if flag
               out (idx+1,cnt,1:idx) = out (idx,i,1:idx);
               out(idx+1,cnt,idx+1) = j;
               cnt = cnt+1;
           end
        end
    end

    idx = idx+1;
    out = compound(out,n,num,d,idx);
end
