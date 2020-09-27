program NACA_5_digit

implicit none
integer::dig1,dig2,dig3,dig45
real::c
!--------------------------------------------------------------------------------------!
write(*,*)"|---------------------------------------------------------------------------|"
write(*,*)"|***************************************************************************|"
write(*,*)"|   Welcome to Foil-Gen-5!                                                  |"
write(*,*)"|   This program generates NACA 5 digit airfoil coordinates                 |" 
write(*,*)"|   Created by Guruprasad S                                                 |"
write(*,*)"|   Email:- guru112358@gmail.com                                            |"
write(*,*)"|   Input the airfoil 5 digit number  as per prompts on the screen          |"
write(*,*)"|   NOTE:- The graphs are plotted using gnuplot                             |"
write(*,*)"|***************************************************************************|"
write(*,*)"|---------------------------------------------------------------------------|"
write(*,*)
write(*,*)
write(*,*)
write(*,*)"Enter the first digit of the airfoil: "
write(*,*)
read *,dig1
write(*,*)
write(*,*)"Enter the second digit of the airfoil: "
read *,dig2
write(*,*)
write(*,*)"Enter the third digit of the airfoil(0 or 1): "
read *,dig3
write(*,*)
write(*,*)"Enter the  last two digits of the airfoil: "
write(*,*)
read *,dig45
write(*,*)
write(*,*)"Enter the desired chord:"
write(*,*)
read *,c
write(*,*)
write(*,*)
!digit 1 -This digit controls the camber. It indicates the designed coefficient of lift (Cl) multiplied by 3/20. In the examble L=2 so Cl=0.3
!digit 2 -The position of maximum camber divided by 20. In the examble P=3 so maximum camber is at 0.15 or 15% chord
!digit 3 - 0 = normal camber line, 1 = reflex camber line
!digits 4 and 5 -The maximum thickness as percentage.In the examble XX=12 so the maximum thickness is 0.12 or 12% chord
!The computation of coordinates is split into separate subroutines for reflexed  and non reflexed airfoils
!the coeffcients scale with cl
if(dig3==0)then
print *,"The airfoil has no reflex"
call Foil_5_gen_nonreflex(dig1,dig2,dig45,c)
else
print *,"The airfoil has a reflex"
call Foil_5_gen_reflex(dig1,dig2,dig3,dig45,c)
end if

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Foil_5_gen_nonreflex(dig1,dig2,dig45,c)
implicit none
integer::dig1,dig2,dig45,i
real::c,p,t,dx,Cl,Cl0,k1,r
real,dimension(0:200)::x,yc,xu,u,xl,yl,yu,theta,dyc_dx,yt

t=(dig45/100.0)*c
p=(dig2)*5.0   !position of maximum camber
dx=c/200.0
 Cl=dig1*0.15
 Cl0=0.30


if(dig2==1)then
k1=361.400*(Cl/Cl0)
r=0.0580*(Cl/Cl0)
else if(dig2==2)then
k1=51.460*(Cl/Cl0)
r=0.1260*(Cl/Cl0)
else if(dig2==3)then
k1=15.957*(Cl/Cl0)
r=0.2025*(Cl/Cl0)
else if(dig2==4)then
k1=6.643*(Cl/Cl0)
r=0.2900*(Cl/Cl0)
else if(dig2==5)then
k1=3.230*(Cl/Cl0)
r=0.3910*(Cl/Cl0)
end if


open(1,file='nonreflex.dat',status='replace')
open(2,file='nonreflex.plt',status='replace')

do i=0,200
x(i)=i*dx
yt(i)=5.0*t*((0.2969*SQRT(x(i)/c))-(0.1260*x(i)/c)-0.3516*((x(i)/c)**2)+(0.2843*(x(i)/c)**3)-(0.1036*(x(i)/c)**4))
if(((x(i)/c)>=0).AND.((x(i)/c)<r))then
yc(i)=(k1/6)*(((x(i)/c)**3)-(3*r*(x(i)/c)**2)+((r**2)*(3-r)*(x(i)/c)))
dyc_dx(i)=(k1/6)*((3*(x(i)/c)**2)-(6*r*(x(i)/c))+((r**2)*(3-r)))
else if(((x(i)/c)>=r).AND.((x(i)/c)<=1))then
yc(i)=((k1*r**3)/6)*(1-(x(i)/c))
dyc_dx(i)=-(k1*r**3)/6
end if


theta(i)=ATAN(dyc_dx(i))
xu(i)=x(i)-yt(i)*SIN(theta(i))
xl(i)=x(i)+yt(i)*SIN(theta(i))

yu(i)=yc(i)+yt(i)*COS(theta(i))
yl(i)=yc(i)-yt(i)*COS(theta(i))

end do

do i=0,200
write(1,*)xu(i),yu(i)
end do

do i=0,200
write(1,*)xl(i),yl(i)
end do
write(*,*)"The airfoil has been generated,dat file can be found in the directory"

 write(2,*)'set xlabel "X"'
 write(2,*)'set ylabel "Y"'
 write(2,*)'set autoscale xy'
 write(2,*)'set title "non reflex 5 series Airfoil surface plot"'
 write(2,*)'set size ratio 0.1'
 write(2,*)'plot "nonreflex.dat" with line'
 CALL SYSTEM('gnuplot -p nonreflex.plt')

 close(1)
 close(2)
 
END SUBROUTINE Foil_5_gen_nonreflex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE Foil_5_gen_reflex(dig1,dig2,dig3,dig45,c)
real::c,p,t,dx,Cl,Cl0,k1,r,k2k1
real,dimension(0:200)::x,yc,xu,u,xl,yl,yu,theta,dyc_dx,yt
integer::dig1,dig2,dig3,dig45,i



t=(dig45/100.0)*c
p=(dig2)*5.0   !position of maximum camber
dx=c/200.0
 Cl=dig1*0.15
 Cl0=0.30

if(dig2==2)then
k1=51.990*(Cl/Cl0)
r=0.1300*(Cl/Cl0)
k2k1=0.000764*(Cl/Cl0)
else if(dig2==3)then
k1=15.793*(Cl/Cl0)
r=0.2170*(Cl/Cl0)
k2k1=0.00677*(Cl/Cl0)
else if(dig2==4)then
k1=6.520*(Cl/Cl0)
r=0.3180*(Cl/Cl0)
k2k1=0.0303*(Cl/Cl0)
else if(dig2==5)then
k1=3.191*(Cl/Cl0)
r=0.4410*(Cl/Cl0)
k2k1=0.1355*(Cl/Cl0)
end if


open(1,file='reflex.dat',status='replace')
open(2,file='reflex.plt',status='replace')

do i=0,200
x(i)=i*dx
yt(i)=5.0*t*((0.2969*SQRT(x(i)/c))-(0.1260*x(i)/c)-0.3516*((x(i)/c)**2)+(0.2843*(x(i)/c)**3)-(0.1036*(x(i)/c)**4))
if(((x(i)/c)>=0).AND.((x(i)/c)<r))then
yc(i)=(k1/6)*(((x(i)/c)-r)**3 -((k2k1*(1-r)**3)*(x(i)/c))-((r**3)*(x(i)/c))+r**3)
dyc_dx(i)=(k1/6)*((3*((x(i)/c)-r)**2)-(k2k1*(1-r)**3)-r**3)
else if(((x(i)/c)>=r).AND.((x(i)/c)<=1))then
yc(i)=(k1/6)*((k2k1*((x(i)/c)-r)**3)-((k2k1*(1-r)**3)*(x(i)/c))-((r**3)*(x(i)/c))+r**3)
dyc_dx(i)=(k1/6)*((3*k2k1*((x(i)/c)-r)**2)-(k2k1*(1-r)**3)-r**3)
end if


theta(i)=ATAN(dyc_dx(i))
xu(i)=x(i)-yt(i)*SIN(theta(i))
xl(i)=x(i)+yt(i)*SIN(theta(i))

yu(i)=yc(i)+yt(i)*COS(theta(i))
yl(i)=yc(i)-yt(i)*COS(theta(i))
end do

do i=0,200
write(1,*)xu(i),yu(i)
end do

do i=0,200
write(1,*)xl(i),yl(i)
end do


write(*,*)"The reflex airfoil has been generated,dat file can be found in the directory"

 write(2,*)'set xlabel "X"'
 write(2,*)'set ylabel "Y"'
 write(2,*)'set autoscale xy'
 write(2,*)'set title "reflex 5 series Airfoil surface plot"'
 write(2,*)'set size ratio 0.1'
 write(2,*)'plot "reflex.dat" with line'
 CALL SYSTEM('gnuplot -p reflex.plt')

 close(1)
 close(2)



END SUBROUTINE Foil_5_gen_reflex




end program NACA_5_digit
