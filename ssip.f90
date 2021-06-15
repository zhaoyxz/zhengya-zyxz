subroutine ssip(a,w,s,m,n,l)    !定义五点平滑子程序 
dimension a(m,n),w(m,n)     !定义数组 
if(l==1)then           
do i=2,m-1       
    do j=2,n-1       
w(i,j)=a(i,j)+s*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4*a(i,j))/4.0    !执行正平滑 
enddo 
enddo 
do i=2,m-1   
    do j=2,n-1   
        a(i,j)=w(i,j) 
    enddo 
Enddo 
return 
Else 
    do i=2,m-1  
        do j=2,n-1 
w(i,j)=a(i,j)+s*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4*a(i,j))/4.0   !执行正平滑  
enddo 
enddo 
do i=2,m-1 
    do j=2,n-1 
        a(i,j)=w(i,j) 
    enddo 
enddo 

do i=2,m-1 
    do j=2,n-1  

w(i,j)=a(i,j)+(-s)*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4*a(i,j))/4.0  !执行正逆平滑 
enddo 
enddo 
do i=2,m-1  
    do j=2,n-1   
        a(i,j)=w(i,j)  
    enddo 
enddo 
endif 
return 
end  !结束子程序