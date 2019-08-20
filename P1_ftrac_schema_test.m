%% iteration avec RK4
for i=1:Niter
    for k=2:Natome
        fext1=forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k);
        Vtemp1(k,:,i+1)=V(k,:,i);
        Ptemp1(k,:,i+1)=P(k,:,i)+dt/2*Vtemp1(k,:,i+1)+(dt/2)^2/2*fext1/m;
    end
    for k=2:Natome       
        fext2(k,:,i+1)=forcetot(Ptemp1(k+1,:,i+1),Ptemp1(k,:,i+1),Ptemp1(:,:,i+1),k)+forcetot(Ptemp1(k-1,:,i+1),Ptemp1(k,:,i+1),Ptemp1(:,:,i+1),k);
        Vtemp2(k,:,i+1)=V(k,:,i)+dt/2*fext2/m;
        Ptemp2(k,:,i+1)=P(k,:,i)+dt/2*Vtemp2(k,:,i+1)+(dt/2)^2/2*fext2/m;
    end
    for k=2:Natome
        fext3(k,:,i+1)=forcetot(Ptemp2(k+1,:,i+1),Ptemp2(k,:,i+1),Ptemp2(:,:,i+1),k)+forcetot(Ptemp2(k-1,:,i+1),Ptemp2(k,:,i+1),Ptemp2(:,:,i+1),k);
        Vtemp3(k,:,i+1)=V(k,:,i)+dt/2*fext3/m;
        Ptemp3(k,:,i+1)=P(k,:,i)+dt*Vtemp3(k,:,i+1)+dt^2/2*fext3/m;
    end
    for k=2:Natome
        fext4(k,:,i+1)=forcetot(Ptemp3(k+1,:,i+1),Ptemp3(k,:,i+1),Ptemp3(:,:,i+1),k)+forcetot(Ptemp3(k-1,:,i+1),Ptemp3(k,:,i+1),Ptemp3(:,:,i+1),k);
        Vtemp4=V(k,:,i)+dt*fext4/m;
    end
    for k=2:Natome
        P(k,:,i+1)=P(k,:,i)+dt/6*(Vtemp1(k,:,i+1)+2*Vtemp2(k,:,i+1)+2*Vtemp4(k,:,i+1)+Vtemp4(k,:,i+1));
    end
    
    fext1_1=forcetot([0,0,0],P(1,:,i),P(:,:,i),1)+forcetot(P(2,:,i),P(1,:,i),P(:,:,i),1);
    Vtemp1_1=V(1,:,i);
    Ptemp1_1=P(1,:,i)+dt/2*Vtemp1_1+(dt/2)^2/2*fext1_1/m;
    fext2_1=forcetot([0,0,0],Ptemp1_1,Ptemp1(:,:,i+1),1)+forcetot(Ptemp1(2,:,i),Ptemp1_1,Ptemp1(:,:,i+1),1);
    Vtemp2_1=V(1,:,i)+dt/2*fext2_1/m;
    Ptemp2_1=P(1,:,i)+dt/2*Vtemp2_1+(dt/2)^2/2*fext2_1/m;
    fext3_1=forcetot([0,0,0],Ptemp2_1,Ptemp2(:,:,i+1),1)+forcetot(Ptemp2(2,:,i),Ptemp2_1,Ptemp2(:,:,i+1),1);
    Vtemp3_1=V(1,:,i)+dt/2*fext3_1/m;
    Ptemp3_1=P(1,:,i)+dt/2*Vtemp3_1+(dt/2)^2/2*fext3_1/m;
    fext4_1=forcetot([0,0,0],Ptemp3_1,Ptemp3(:,:,i+1),1)+forcetot(Ptemp3(2,:,i),Ptemp3_1,Ptemp3(:,:,i+1),1);
    Vtemp4_1=V(1,:,i)+dt*fext4_1/m;
    P(1,:,i+1)=P(1,:,i)+dt/6*(Vtemp1_1+2*Vtemp2_1+2*Vtemp3_1+Vtemp4_1);
    
    fext1_d=forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1)+ftrac;
    Vtemp1_d=V(Natome+1,:,i);
    Ptemp1_d=P(Natome+1,:,i)+dt/2*Vtemp1_d+(dt/2)^2/2*fext1_d/m;
    fext2_d=forcetot(Ptemp1(Natome,:,i),Ptemp1_d,Ptemp1(:,:,i),Natome+1)+ftrac;
    Vtemp2_d=V(Natome+1,:,i)+dt/2*fext2_d/m;
    Ptemp2_d=P(Natome+1,:,i)+dt/2*Vtemp2_d+(dt/2)^2/2*fext2_d/m;
    fext3_d=forcetot(Ptemp2(Natome,:,i),Ptemp2_d,Ptemp2(:,:,i),Natome+1)+ftrac;
    Vtemp3_d=V(Natome+1,:,i)+dt/2*fext3_d/m;
    Ptemp3_d=P(Natome+1,:,i)+dt/2*Vtemp3_d+(dt/2)^2/2*fext3_d/m;
    fext4_d=forcetot(Ptemp3(Natome,:,i),Ptemp3_d,Ptemp3(:,:,i),Natome+1)+ftrac;
    Vtemp4_d=V(Natome+1,:,i)+dt*fext4_d/m;
    P(Natome+1,:,i+1)=P(Natome+1,:,i)+dt/6*(Vtemp1_d+2*Vtemp2_d+2*Vtemp3_d+Vtemp4_d);
end
    


