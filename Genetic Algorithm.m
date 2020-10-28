clear
format short;

upper = 10.0;
lower = -10.0;
num = 1000;
cross_P= 0.99;
precision = 6;
% 取亂數
for i = 1:num
    for j= 1:4
        x(i,j) = -10 + (10-(-10))*rand;
    end
end
max_cycle = 500;


x(:,:) = fix(x(:,:)*(10^precision))/(10^precision);  %取小數點以下6位

for cycle = 1:max_cycle
    % fx
    x(:,5) = 100.*(x(:,2)-x(:,1).^2).^2 + (1-x(:,1)).^2 + 90.*(x(:,4)-x(:,3).^2).^2+...
        (1-x(:,3)).^2 + 10.1.*((x(:,2)-1).^2 + (x(:,4)-1).^2) + 19.8.*(x(:,2)-1).*(x(:,4)-1);
    
    % 整理為只小數點後6位
    x(:,:) = fix(x(:,:)*(10^precision))/(10^precision);
    
    min_F(cycle,5) = min(x(:,5));
    z = find(x(:,5)==min(min(x(:,5))));  %找出index
    min_F(cycle,1) = x(z(1,1),1);
    min_F(cycle,2) = x(z(1,1),2);          %index 的第一個
    min_F(cycle,3) = x(z(1,1),3);
    min_F(cycle,4) = x(z(1,1),4);
    fprintf(' cycle%3d    X1 = %-3.6f  X2 = %-3.6f  X3 = %-3.6f  X4 = %-3.6f  fitness = %-3.6f\n',cycle,min_F(cycle,1),min_F(cycle,2),min_F(cycle,3),min_F(cycle,4),min_F(cycle,5));
    
    
    %排序後暫存
    [B,index] = sortrows(x(:,5));
    for i = 1:num
        for j = 1:5
            temp(i,j) = x(index(i),j);
        end
    end
    
    %暫存帶回
    x(:,:) = temp(:,:);
    clear temp;
    
    %被選機率+累進機率
    x(:,6) = x(:,5)./sum(x(:,5));
    for i = 1:num
        if i > 1
            x(i,7) = x(num+1-i,6)+x(i-1,7);
        else
            x(i,7) = x(num+1-i,6);
        end
    end
    
    %取亂數
    for i =1:num
        random(i,1) = rand;
     end
    
    %找出新群體位置
    for i= 1:num
        for j= 1:num
            if random(i,1) < x(j,7)
                random(i,2) = j;
                break;
            end
        end
    end
    
    
    %輪盤染色體
    for i=1:num
        for j=1:4
            temp(i,j) = x(random(i,2),j);
        end
    end
    
    %交叉亂數
    for i=1:num
        random(i,3) = rand;
    end
    
    %找需要交叉的基因列
    k = 0;
    for i= 1:num
        if random(i,3)< cross_P
            k=k+1;
            random(k,4) = i;
        end
    end
    
    %交叉
    for i = 1:k
        if mod(i,2) ==0
            cross_R = fix(4*rand+1);  % 取1 ~ n 的亂數為交叉點   fix為整理為整數
            beta = rand;
            for j= 1:cross_R-1
                temp(num+i-1,j) = x(random(random(i-1,4),2),j);
                temp(num+i,j) = x(random(random(i,4),2),j);
            end
            Xk = temp(random(random(i-1,4),2),cross_R);
            Yk = temp(random(random(i,4),2),cross_R);
            temp(num+i-1,cross_R) = beta*Xk+(1-beta)*Yk;
%             temp(num+i,cross_R) = beta*Yk+(1-beta)*Xk;
            temp(num+i,cross_R) = lower + beta*(upper-lower);
            for j = cross_R+1 : 4
                temp(num+i-1,j) = x(random(random(i,4),2),j);  
                temp(num+i,j) =  x(random(random(i-1,4),2),j);
            end
        end
    end
    
    %暫存基因檔
    mutation(:,:) = temp(:,:);
    
    %增加num的總數方便後面計算
    if mod(k,2)==1
        new_num = num+ k-1;
    else
        new_num = num+k;
    end
    
    %突變 + 新增辨識列，辨識是否有更新
    pm = 0.01;
    for i=1:new_num
        test = 0;
        for j= 1:4
            mutation_R = rand;
            if mutation_R < pm
                mutation_R = rand;
                beta  =rand;
                if mutation_R< 0.5
                    y = upper-mutation(i,j);
                    mutation(i,j) = mutation(i,j)+y*beta*(1-(cycle/max_cycle));
                    
                else
                    y = mutation(i,j)-lower;
                    mutation(i,j) = mutation(i,j)-y*beta*(1-(cycle/max_cycle));
                end
                test = test+1;
            end
        end
        temp(i,5) = test;
    end
    
    %將有突變的原基因放到最後面
    x = 1;
    for i = 1:new_num
        if temp(i,5)>0
            for j = 1:4
                mutation(new_num+x,j) = temp(i,j);
            end
            x = x +1;
        end
    end
    new_num = new_num + x-1;
    
    % fx
    mutation(:,5) = 100.*(mutation(:,2)-mutation(:,1).^2).^2 + (1-mutation(:,1)).^2 + 90.*(mutation(:,4)-mutation(:,3).^2).^2+...
        (1-mutation(:,3)).^2 + 10.1.*((mutation(:,2)-1).^2 + (mutation(:,4)-1).^2) + 19.8.*(mutation(:,2)-1).*(mutation(:,4)-1);
    
    if new_num >= num
        [B,index] = sortrows(mutation(:,5));
        new_num = num;
    end
    
    for i=1: num*0.01
        for j = 1:4
            x(i,j) = mutation(index(i),j);
        end
    end
    
    for i = 1:num*0.99
        for j= 1:4
            x(num*0.01+i,j) = -10 + (10-(-10))*rand;
        end
    end
    x(:,:) = fix(x(:,:)*(10^precision))/(10^precision);
    
    %清空
    clear mutation;
    clear random;
    clear B;
    clear index;
    clear temp;
    
end 

% figure('Name','value');
subplot(5,1,1),plot(min_F(:,5),'linewidth',1.5);
title('fitness');
subplot(5,1,2),plot(min_F(:,1),'linewidth',1.5);
title('X1');
subplot(5,1,3),plot(min_F(:,2),'linewidth',1.5);
title('X2');
subplot(5,1,4),plot(min_F(:,3),'linewidth',1.5);
title('X3');
subplot(5,1,5),plot(min_F(:,4),'linewidth',1.5);
title('X4');









