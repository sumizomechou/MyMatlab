    H=[155,160,165];        %�ŵ�����
    S=zeros(1,903);        %״̬
    ACTION=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];   %������
    gamma=0.95;
    Ts=1;   %ÿ��ʱ�ڵĴ���ʱ��
    symbol=0;   %�������Ԫ��
    allsymbol=0; 
    epsilon=0.1;      %̽����
    alpha=0.1;
    eta=0.2;
    q_value=zeros(903,15);   %��ʼ��Q��
    q1_value=zeros(903,15);
 for i=1:50000
    E=600;   %�����������
    E1=600;  %�������ʣ������
    b=randperm(40);  %��ʼ����ʼ״̬
    c=b(1);  %ȡb�ĵ�һ��
    if c<=2
        h=155;  %�ŵ�����
    elseif c>=3 && c<=21
        h=160;
    else
        h=165;
    end
    switch h
        case 155
            state=301;   %301-1
        case 160
            state=602;   %602-302
        case 165
            state=903;   %903-603
    end
    while E1>0
       if binornd(1,epsilon)==1     %epsilon-̰����ѡ����
           b=ACTION(randperm(length(ACTION)));
           action=b(1); 
       else
           b=q_value(state,:);
           [~,n]=max(b);
           action=n;
       end
       switch action  %��ͬ���ƶ�Ӧ�Ĺ���
           case {1,6,11}  % 1 2345 6 78910 11 12131415
               P=2;
           case {2,7,12}  %1 2 345 6 7 8910 11 12 131415
               P=4;
           case {3,8,13}
               P=6;
           case {4,9,14}
               P=8;
           case {5,10,15}  %15=3�ֵ���x5������
               P=10;
       end
        if (P*Ts>E1)
            switch E1       %�ڿ�ѡ�Ķ�����ʹ��epsilon-̰����ѡ����
                case 8
                    ACTION=[1 2 3 4 5 6 7 8 9 10 11 12];
                    if binornd(1,epsilon)==1     
                       b=ACTION(randperm(length(ACTION)));
                       n=b(1);
                       action=n;
                    else
                        q1=q_value(:,1:4);
                        q2=q_value(:,6:9);
                        q3=q_value(:,11:14);
                        q=[q1,q2,q3];
                        b=q(state,:);
                        [~,n]=max(b);
                        action=n;
                    end
                    switch action
                        case {1,5,9}
                            P=2;
                        case {2,6,10}
                            P=4;
                        case {3,7,11}
                            P=6;
                        case {4,8,12}
                            P=8;
                    end
                    if n<5
                        action=n;
                    elseif n>4 && n<9
                        action=n+1;
                    else
                        action=n+2;
                    end
                case 6
                    ACTION=[1 2 3 4 5 6 7 8 9];
                    if binornd(1,epsilon)==1     
                       b=ACTION(randperm(length(ACTION)));
                       n=b(1); 
                       action=n;
                    else
                        q1=q_value(:,1:3);
                        q2=q_value(:,6:8);
                        q3=q_value(:,11:13);
                        q=[q1,q2,q3];
                        b=q(state,:);
                        [~,n]=max(b);
                        action=n;
                    end
                    switch action
                        case {1,4,7}
                            P=2;
                        case {2,5,8}
                            P=4;
                        case {3,6,9}
                            P=6;
                    end
                    if n<4
                        action=n;
                    elseif n>3 && n<7
                        action=n+2;
                    else
                        action=n+4;
                    end
                case 4
                    ACTION=[1 2 3 4 5 6];
                    if binornd(1,epsilon)==1     
                       b=ACTION(randperm(length(ACTION)));
                       n=b(1); 
                       action=n;
                    else
                        q1=q_value(:,1:2);
                        q2=q_value(:,6:7);
                        q3=q_value(:,11:12);
                        q=[q1,q2,q3];
                        b=q(state,:);
                        [~,n]=max(b);
                        action=n;
                    end
                    switch action
                        case {1,3,5}
                            P=2;
                        case {2,4,6}
                            P=4;
                    end
                    if n<3
                        action=n;
                    elseif n>2 && n<5
                        action=n+3;
                    else
                        action=n+6;
                    end
                case 2
                    ACTION=[1 2 3];
                    if binornd(1,epsilon)==1     
                       b=ACTION(randperm(length(ACTION)));
                       n=b(1);
                    else
                        q1=q_value(:,1);
                        q2=q_value(:,6);
                        q3=q_value(:,11);
                        q=[q1,q2,q3];
                        b=q(state,:);
                        [~,n]=max(b);
                    end
                    P=2;
                    if n==1
                        action=n;
                    elseif n==2
                        action=n+4;
                    else
                        action=n+8;
                    end
                case 0
                    P=0;
           end
        end

       SL=171+10*log10(P)+10*log10(eta);  %165��175֮��
       snr=SL-h;
      switch action
          case {1,2,3,4,5}
              BER=0.5*erfc(sqrt(snr));   %������
              M=2;   %���ƽ���
          case {6,7,8,9,10}
              BER=erfc(sqrt(4*snr)*sin(pi/8));
              M=4;
          case {11,12,13,14,15}
              BER=erfc(sqrt(6*snr)*sin(pi/16));
              M=8;
      end

          R=10*(0.05-BER)*M/P;
      
      symbol=symbol+M*(1-BER);
      allsymbol=allsymbol+M;
      E1=E1-P*Ts;

     c=randperm(40);    %�����µ��ŵ�
     d=c(1);
     if d<=2
        next_h=155;
     elseif d>=3 && d<=21
        next_h=160;
     else
        next_h=165;
     end
      switch next_h         %�����ŵ�״̬����ʣ����
        case 155
                next_state=E1/2+1;   %598/2+302=300
        case 160
                next_state=E1/2+302;   %598/2+302=601
        case 165
                next_state=E1/2+603;   %598/2+603=902
      end
      if E1>8           %�ڿ�ѡ�Ķ�������̰����ѡ������ֵ�Ķ���
         d=q_value(next_state,:);
         [~,n]=max(d);
         max_action=n;
      else
         switch E1
             case 8
                 q1=q_value(:,1:4);
                 q2=q_value(:,6:9);
                 q3=q_value(:,11:14);
                 q=[q1,q2,q3];
                 d=q(next_state,:);
                 [~,n]=max(d);
                 if n<5
                        max_action=n;
                    elseif n>4 && n<9
                        max_action=n+1;
                    else
                        max_action=n+2;
                 end 
             case 6
                 q1=q_value(:,1:3);
                 q2=q_value(:,6:8);
                 q3=q_value(:,11:13);
                 q=[q1,q2,q3];
                 d=q(next_state,:);
                 [~,n]=max(d);
                 if n<4
                        max_action=n;
                    elseif n>3 && n<7
                        max_action=n+2;
                    else
                        max_action=n+4;
                 end
             case 4
                 q1=q_value(:,1:2);
                 q2=q_value(:,6:7);
                 q3=q_value(:,11:12);
                 q=[q1,q2,q3];
                 d=q(next_state,:);
                 [~,n]=max(d);
                 if n<3
                        max_action=n;
                    elseif n>2 && n<5
                        max_action=n+3;
                    else
                        max_action=n+6;
                  end
             case 2
                 q1=q_value(:,1);
                 q2=q_value(:,6);
                 q3=q_value(:,11);
                 q=[q1,q2,q3];
                 d=q(next_state,:);
                 [~,n]=max(d);
                 if n==1
                        max_action=n;
                    elseif n==2
                        max_action=n+4;
                    else
                        max_action=n+8;
                 end
         end
      end
      q_value(state,action)=q_value(state,action)+alpha*(R+gamma* q_value(next_state,max_action)-q_value(state,action));   
      state=next_state;
      h=next_h;
    end
 end
  wumalv=1-symbol/allsymbol;
  Xa = sprintf('�޸ĺ�:\n������Ԫ%d�����ͳɹ���Ԫ%d',allsymbol,symbol);
  Xb = sprintf('������%d',wumalv);
  disp(Xa)
  disp(Xb)
 %save q_table_channel_1.mat q_value   %�������ɵ�Q��