!%%%%%%%%%%%%%%%
!Navier-Stokes's equation solver in spherical coordinates with Boussinesq approx
!To ommit convection effects, set gravity values to zero, thereby external forces are null
!%%%%%%%%%%%%%%%
program NS3DTemp
implicit none

integer, parameter :: nr=30,nth=20,nphi=30

integer :: i,j,k,t,Tmax,pt,Nplot,NmaxP,&
	imin,imax,jmin,jmax,kmin,kmax,cont
real :: lth,lr,lphi,Pi,Rmin,Rmax,dth,dr,dphi,Re,Ga,dt,&
	ti,tf2,beta,Pr,gr,gth,gphi

real, dimension(:),allocatable :: th,r,phi
real, dimension(:,:),allocatable :: T0
real, dimension(:,:,:),allocatable:: &
	u1,u2,v1,v2,w1,w2,p1,p2,F,G,L,flth,flr,flphi,TE1,TE2

imin=1 ; imax=nth
jmin=1 ; jmax=nr
kmin=1 ; kmax=nphi

ALLOCATE(th(imin-1:imax+1),r(jmin-1:jmax+1),&
	phi(kmin-1:kmax+1))
ALLOCATE(U1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	U2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	V1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	V2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	W1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	W2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	P1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	P2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	F(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	G(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	L(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
flth(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
flr(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
flphi(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(TE1(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),&
	TE2(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1))
ALLOCATE(T0(imin-1:imax+1,kmin-1:kmax+1))

Re=4560.0
Ga=0.9

Pr=0.71
beta=3.625e-3
gr=-1.2
gth=0.0
gphi=0.0

Pi=ACOS(-1.0)
Nplot=500

Rmin=0.5 ; Rmax=1.0
lr=Rmax-Rmin ; lth=2.0*Pi ; lphi=Pi

dr=lr/float(nr)
dth=lth/float(nth)
dphi=lphi/float(nphi)


!Time parameters
dt=0.005*(Re*Pr/2.0)/( (1.0/(dr*dr)) + (1.0/(dphi*dphi)) + (1.0/(dth*dth)) )
Tmax=500000
print*, 'dt=',dt,' Tmax=', Tmax

!---------------- GRID --------------------------------
Th(imin-1)=0.
Th(imin)=0.5*dth
DO i=(imin+1),imax
Th(i)=0.5*dth+(i-imin)*dth
END DO
Th(imax+1)=Lth

r(jmin-1)=Rmin
r(jmin)=Rmin+0.5*dr
DO J=(jmin+1),jmax
r(J)=Rmin+0.5*dr+(j-jmin)*dr
END DO
r(jmax+1)=Rmax

Phi(kmin-1)=0.0
Phi(kmin)=0.5*dphi
DO k=(kmin+1),kmax
Phi(k)=0.5*dphi+(k-kmin)*dphi
END DO
Phi(kmax+1)=Lphi
!-------------------------------------------------------

!Temp at inner sphere
DO i=imin-1,imax+1
DO k=kmin-1,kmax+1
T0(i,k)=1.0!*sin(phi(k))
END DO
END DO

u1=0.0 ; v1=0.0 ; w1=0.0 ; p1=0.0 ; TE1=0.0
u2=u1 ;v2=v1 ; w2=w1 ; p2=p1 ; TE2=TE1

cont=0
call CPU_time(ti)
DO t=1, Tmax
CALL BOUNDARY(u1,u2,v1,v2,w1,w2,F,G,L,TE1,TE2,T0,&
	imin,imax,jmin,jmax,kmin,kmax)

CALL TEMPERATURE(TE1,TE2,u1,v1,w1,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,Pr,dr,dth,dphi,dt)

CALL FCALC(F,TE2,u1,v1,w1,flth,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,dr,dth,dphi,dt,gth,beta)

CALL GCALC(G,TE2,u1,v1,w1,flr,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,dr,dth,dphi,dt,gr,beta)

CALL LCALC(L,TE2,u1,v1,w1,flphi,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,dr,dth,dphi,dt,gphi,beta)

CALL PRESSURE(P1,P2,F,G,L,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	dth,dr,dphi,dt)


CALL U2V2W2(u2,v2,w2,F,G,L,p2,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	dth,dphi,dr,dt)

IF(mod(t,Nplot).EQ.0) THEN

IF(cont.eq.0) then
call CPU_time(tf2)
write(*,*) 'Tiempo estimado=',((Tmax-Nplot)/Nplot)*(tf2-ti)/60.0, ' min'
cont=1
end if

 1    format(1x,'N=',I8,'  U=',F12.4,'  V=',F12.4,' W=',F12.4,' P=',F12.4, ' T=', F12.4 )

	print 1,t,U2(nth/2,nr/2,nphi/2),V2(nth/2,nr/2,nphi/2),W2(nth/2,nr/2,nphi/2),P2(nth/2,nr/2,nphi/2), TE2(nth/2,nr/2,nphi/2)

!CALL ANIMACION(r,th,phi,u2,v2,w2,&
!	imin,imax,jmin,jmax,kmin,kmax,Nplot,t,TE2)
END IF

u1=u2 ; v1=v2 ; w1=w2 ; TE1=TE2
ENDDO

do k=kmin,kmax
do j=jmin,jmax
do i=imin,imax
u2(i,j,k)=0.5*(u2(i,j,k)+u2(i-1,j,k))
v2(i,j,k)=0.5*(v2(i,j,k)+v2(i,j-1,k))
w2(i,j,k)=0.5*(w2(i,j,k)+w2(i,j,k-1))
enddo
enddo
enddo

2    format(2x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5)
OPEN(1,file='data.dat')

write(1,*) ' TITLE = "NS2DSphere" '
write(1,*) ' VARIABLES = "X", "Y","Z", "U", "V", "W","MAGU","TEMP"'
write(1,*) ' ZONE T="1", I=',imax+2,', J=',jmax+2,', K=',kmax

DO k=kmin,kmax
DO j=jmin-1,jmax+1
DO i=imin-1,imax+1

!DO i=imin-1,imax+1 !For paraview

	write(1,2) r(j)*Sin(phi(k))*Cos(th(i)), &
		r(j)*Sin(phi(k))*Sin(th(i)), &
		r(j)*Cos(phi(k)), &
		V2(I,J,k)*Cos(Th(i))*Sin(phi(k))-U2(I,J,k)*Sin(Th(i))+W2(I,J,k)*Cos(Th(i))*Cos(phi(k)), &
		V2(I,J,k)*Sin(Th(i))*Sin(phi(k))+U2(I,J,k)*Cos(Th(i))+W2(I,J,k)*Sin(Th(i))*Cos(phi(k)), &
		V2(I,J,k)*Cos(phi(k))-W2(I,J,k)*Sin(phi(k)),&
		SQRT((V2(I,J,k)*Cos(Th(i))*Sin(phi(k))-U2(I,J,k)*Sin(Th(i))+W2(I,J,k)*Cos(Th(i))*Cos(phi(k)))**2&
		+(V2(I,J,k)*Sin(Th(i))*Sin(phi(k))+U2(I,J,k)*Cos(Th(i))+W2(I,J,k)*Sin(Th(i))*Cos(phi(k)))**2&
		+(V2(I,J,k)*Cos(phi(k))-W2(I,J,k)*Sin(phi(k)))**2),TE2(i,j,k)
ENDDO
ENDDO
ENDDO
CLOSE(1)


END PROGRAM NS3DTemp

!%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%

subroutine BOUNDARY(u1,u2,v1,v2,w1,w2,F,G,L,TE1,TE2,T0,&
	imin,imax,jmin,jmax,kmin,kmax)
implicit none
integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, dimension(imin-1:imax+1,kmin-1:kmax+1), intent(inout) :: T0
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: u1,u2,v1,v2,w1,w2,F,G,L,TE1,TE2
integer :: i,j,k

!j=jmin-1 ------ Rmin
FORALL(i=imin:imax,k=kmin:kmax)
v1(i,jmin-1,k)=0.0
v2(i,jmin-1,k)=0.0
u1(i,jmin-1,k)=2.*(0.0)-u1(i,jmin,k)
u2(i,jmin-1,k)=2.*(0.0)-u2(i,jmin,k)
w1(i,jmin-1,k)=-w1(i,jmin,k)
w2(i,jmin-1,k)=-w2(i,jmin,k)
G(i,jmin-1,k)=v1(i,jmin-1,k)
TE1(i,jmin-1,k)=2.*(T0(i,k))-TE1(i,jmin,k)
TE2(i,jmin-1,k)=2.*(T0(i,k))-TE2(i,jmin,k)
END FORALL

!j=jmax+1 ------- Rmax
FORALL(i=imin:imax,k=kmin:kmax)
v1(i,jmax,k)=0.
v2(i,jmax,k)=0.
u1(i,jmax+1,k)=-u1(i,jmax,k)
u2(i,jmax+1,k)=-u2(i,jmax,k)
w1(i,jmax+1,k)=-W1(i,jmax,k)
w2(i,jmax+1,k)=-W2(i,jmax,k)
G(i,jmax,k)=v2(i,jmax,k)
TE1(i,jmax+1,k)=2.*(0.0)-TE1(i,jmax,k)
TE2(i,jmax+1,k)=2.*(0.0)-TE2(i,jmax,k)
END FORALL

!------- PERIODIC ----------
FORALL(j=jmin:jmax,k=kmin:kmax)
v1(imin-1,j,k)=v1(imax-1,j,k)
v2(imin-1,j,k)=v2(imax-1,j,k)
u1(imin-1,j,k)=u1(imax-1,j,k)
u2(imin-1,j,k)=u2(imax-1,j,k)
w1(imin-1,j,k)=w1(imax-1,j,k)
w2(imin-1,j,k)=w2(imax-1,j,k)

v1(imax,j,k)=v1(imin,j,k)
v2(imax,j,k)=v2(imin,j,k)
u1(imax,j,k)=u1(imin,j,k)
u2(imax,j,k)=u2(imin,j,k)
w1(imax,j,k)=w1(imin,j,k)
w2(imax,j,k)=w2(imin,j,k)

v1(imax+1,j,k)=v1(imin+1,j,k)
v2(imax+1,j,k)=v2(imin+1,j,k)
u1(imax+1,j,k)=u1(imin+1,j,k)
u2(imax+1,j,k)=u2(imin+1,j,k)
w1(imax+1,j,k)=w1(imin+1,j,k)
w2(imax+1,j,k)=w2(imin+1,j,k)

F(imin-1,j,k)=u2(imin-1,j,k)
F(imax,j,k)=u2(imax,j,k)

TE1(imin-1,J,k)=TE1(imin,J,k)
TE2(imin-1,J,k)=TE2(imin,J,k)
TE1(imax+1,J,k)=TE1(imax,J,k)
TE2(imax+1,J,k)=TE2(imax,J,k)
TE1(imin,J,k)=TE1(imax,J,k)
TE2(imin,J,k)=TE2(imax,J,k)
END FORALL

!k=kmax+1 -----------
FORALL(i=imin:imax,j=jmin:jmax)
u1(i,j,kmax+1)=u1(i,j,kmax)
u2(i,j,kmax+1)=u2(i,j,kmax)
v1(i,j,kmax+1)=v1(i,j,kmax)
v2(i,j,kmax+1)=v2(i,j,kmax)
w1(i,j,kmax)=w1(i,j,kmax-1)
w2(i,j,kmax)=w2(i,j,kmax-1)
L(i,j,kmax)=w2(i,j,kmax)
TE1(I,J,kmax+1)=TE1(I,J,kmax)
TE2(I,J,kmax+1)=TE2(I,J,kmax)
END FORALL

!k=kmin-1 ------------
FORALL(i=imin:imax,j=jmin:jmax)
u1(i,j,kmin-1)=u1(i,j,kmin)
u2(i,j,kmin-1)=u2(i,j,kmin)
v1(i,j,kmin-1)=v1(i,j,kmin)
v2(i,j,kmin-1)=v2(i,j,kmin)
w1(i,j,kmin-1)=w1(i,j,kmin)
w2(i,j,kmin-1)=w2(i,j,kmin)
L(i,j,kmin-1)=W2(i,j,kmin-1)
TE1(I,J,kmin-1)=TE1(I,J,kmin)
TE2(I,J,kmin-1)=TE2(I,J,kmin)
END FORALL

return
end subroutine BOUNDARY

subroutine TEMPERATURE(TE1,TE2,u1,v1,w1,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,Pr,dr,dth,dphi,dt)
implicit none

integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: Re,Ga,Pr,dr,dth,dphi,dt
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: TE1,TE2,u1,v1,w1
real, intent(inout) :: r(jmin-1:jmax+1),&
			phi(kmin-1:kmax+1)

integer :: i,j,k
real :: d2Tdr2,dTdr,d2Tdphi2,dTdphi,&
	d2Tdth2,dTvdr,dTudth,dTwdphi

do i=imin,imax
do j=jmin,jmax
do k=kmin,kmax
d2Tdr2=(TE1(i,j+1,k)-2.*TE1(i,j,k)+TE1(i,j-1,k))/(dr*dr)
dTdr=(TE1(i,j+1,k)-TE1(i,j-1,k))/(2.*dr)
d2Tdphi2=(TE1(i,j,k+1)-2.*TE1(i,j,k)+TE1(i,j,k-1))/(dphi*dphi)
dTdphi=(TE1(i,j,k+1)-TE1(i,j,k-1))/(2.*dphi)
d2Tdth2=(TE1(i+1,j,k)-2.*TE1(i,j,k)+TE1(i-1,j,k))/(dth*dth)

DTUDth= &
(1./dth)* &
(U1(I,J,k)*((TE1(I,J,k)+TE1(I+1,J,k))/2.) - &
U1(I-1,J,k)*((TE1(I-1,J,k)+TE1(I,J,k))/2.))+ &
(Ga/dth)* &
(ABS(U1(I,J,k))*((TE1(I,J,K)-TE1(I+1,J,K))/2.)- &
ABS(U1(I-1,J,K))*((TE1(I-1,J,K)-TE1(I,J,K))/2.))

DTVDr= &
(1./dr)* &
(V1(I,J,K)*((TE1(I,J,K)+TE1(I,J+1,K))/2.)- &
V1(I,J-1,K)*((TE1(I,J-1,K)+TE1(I,J,K))/2.))+ &
(Ga/dr)* &
(ABS(V1(I,J,K))*((TE1(I,J,K)-TE1(I,J+1,K))/2.)- &
ABS(V1(I,J-1,K))*((TE1(I,J-1,K)-TE1(I,J,K))/2.))

DTWDphi= &
(1./dphi)* &
(W1(I,J,K)*((TE1(I,J,K)+TE1(I,J,K+1))/2.) - &
W1(I,J,K-1)*((TE1(I,J,K-1)+TE1(I,J,K))/2.))+ &
(Ga/dphi)* &
(ABS(W1(I,J,K))*((TE1(I,J,K)-TE1(I,J,K+1))/2.)- &
ABS(W1(I,J,K-1))*((TE1(I,J,K-1)-TE1(I,J,K))/2.))

TE2(i,j,k)=TE1(i,j,k)+dt*(&
	(1./(Pr*Re))*(d2Tdr2+(1./r(j))*dTdr + (1./(r(j)*r(j)))*d2Tdphi2 +&
	dTdphi*(1./(tan(phi(K)))) +d2Tdth2*(1./(r(j)*r(j)*sin(phi(K))*sin(phi(K)))))&
	-dTVdr-dTUdth*(1./(r(j)*sin(phi(K))))-dTwdphi*(1./r(j))-2.*TE1(i,j,k)*V1(i,j,k)/r(j)&
	- TE1(i,j,k)*W1(i,j,k)*(1./(r(j)*tan(phi(k)))) )
enddo
enddo
enddo

return
end subroutine TEMPERATURE

subroutine FCALC(F,TE2,u1,v1,w1,flth,r,phi,&
		imin,imax,jmin,jmax,kmin,kmax,&
		Re,Ga,dr,dth,dphi,dt,gth,beta)
implicit none

integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: Re,Ga,dr,dth,dphi,dt,gth,beta
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: F,TE2,u1,v1,w1,flth
real, intent(inout) :: r(jmin-1:jmax+1),&
			phi(kmin-1:kmax+1)

integer :: i,j,k
real :: d2udr2,dudr,d2udphi2,dudphi,d2udth2,&
	dvdth,dwdth,dvudr,du2dth,dwudphi

real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1) :: Lulu

do i=imin,imax-1
do j=jmin,jmax
do k=kmin,kmax
d2udr2=(u1(i,j+1,k)-2.*u1(i,j,k)+u1(i,j-1,k))/(dr*dr)
dudr=(u1(i,j+1,k)-u1(i,j-1,k))/(2.*dr)
d2udphi2=( u1(i,j,k+1)-2.*u1(i,j,k)+u1(i,j,k-1) )/(dphi*dphi)
dudphi=(u1(i,j,k+1)-u1(i,j,k-1))/(2.*dphi)
d2udth2=(u1(i+1,j,k)-2.*u1(i,j,k)+u1(i-1,j,k))/(dth*dth)
dvdth=(v1(i+1,j,k)-v1(i-1,j,k))/(2.*dth)
dwdth=(w1(i+1,j,k)-w1(i-1,j,k))/(2.*dth)
dvudr=0.25/dr*((V1(I,J,K)+V1(I+1,J,K))*(U1(I,J,K)+U1(I,J+1,K))- &
	(V1(I,J-1,K)+V1(I+1,J-1,K))*(U1(I,J-1,K)+U1(I,J,K)))+ &
	0.25*Ga/dr* &
	(ABS(V1(I,J,K)+V1(I+1,J,K))*(U1(I,J,K)-U1(I,J+1,K))- &
	ABS(V1(I,J-1,K)+V1(I+1,J-1,K))*(U1(I,J-1,K)-U1(I,J,K)))
du2dth=0.25/dth*((U1(I,J,K)+U1(I+1,J,K))**2-(U1(I-1,J,K)+U1(I,J,K))**2)+ &
	0.25*Ga/dth* &
	(ABS(U1(I,J,K)+U1(I+1,J,K))*(U1(I,J,K)-U1(I+1,J,K))- &
	ABS(U1(I-1,J,K)+U1(I,J,K))*(U1(I-1,J,K)-U1(I,J,K)))
dwudphi=0.25/dphi*((W1(I,J,K)+W1(I+1,J,K))*(U1(I,J,K)+U1(I,J,K+1))- &
	(W1(I,J,K-1)+W1(I+1,J,K-1))*(U1(I,J,K-1)+U1(I,J,K)))+ &
	0.25*Ga/dphi* &
	(ABS(W1(I,J,K)+W1(I+1,J,K))*(U1(I,J,K)-U1(I,J,K+1))- &
	ABS(W1(I,J,K-1)+W1(I+1,J,K-1))*(U1(I,J,K-1)-U1(I,J,K)))

LULU(I,J,k)= &
	D2UDr2+(2./r(j))*DUDr &
	+(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*D2UDth2 &
	+(1./(r(j)*r(j)))*D2UDphi2 &
	+(1./(r(j)*r(j)*Tan(phi(k))))*DUDphi


flth(i,j,k)=(1. - 0.5*beta*(TE2(i+1,j,k)+TE2(i,j,k)))*gth

F(I,J,k)= &
	U1(I,J,k)+ &
	dt*( &
	(1./Re)*( &
	LULU(I,J,k)+(2./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*DVDth &
	+(2./(r(j)*r(j)))*(Cos(phi(k))/(Sin(phi(k))*Sin(phi(k))))*DWDth &
	-(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*U1(I,J,k) &
	) &        
	-DVUDr-(1./(r(j)*Sin(phi(k))))*DU2Dth &
	-(1./r(j))*DWUDphi &
	-(3./r(j))*(U1(I,J,k)*V1(I,J,k)) &
	-(2./(r(j)*Tan(phi(k))))*(U1(I,J,k)*W1(I,J,k)) &
	+FLth(I,J,k) &
	)
enddo
enddo
enddo

return
end subroutine FCALC

subroutine GCALC(G,TE2,u1,v1,w1,flr,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,dr,dth,dphi,dt,gr,beta)
implicit none
integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: Re,Ga,dr,dth,dphi,dt,gr,beta
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: G,TE2,u1,v1,w1,flr
real, intent(inout) :: r(jmin-1:jmax+1),&
			phi(kmin-1:kmax+1)

integer :: i,j,k
real :: d2vdr2,dvdr,d2vdphi2,dvdphi,d2vdth2,&
	dwdphi,dudth,dv2dr,duvdth,dwvdphi

real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1) :: Lv

do i=imin,imax
do j=jmin,jmax-1
do k=kmin,kmax

d2vdr2=(v1(i,j+1,k)-2.*v1(i,j,k)+v1(i,j-1,k))/(dr*dr)
dvdr=(v1(i,j+1,k)-v1(i,j-1,k))/(2.*dr)
d2vdphi2=(v1(i,j,k+1)-2.*v1(i,j,k)+v1(i,j,k-1) )/(dphi*dphi)
dvdphi=(v1(i,j,k+1)-v1(i,j,k-1))/(2.*dphi)
d2vdth2=(v1(i+1,j,k)-2.*v1(i,j,k)+v1(i-1,j,k) )/(dth*dth)
dwdphi=(w1(i,j,k+1) - w1(i,j,k-1))/(2.*dphi)
dudth=(u1(i+1,j,k)-u1(i-1,j,k))/(2.*dth)
dv2dr=0.25/dr*((V1(I,J,K)+V1(I,J+1,K))**2-(V1(I,J-1,K)+V1(I,J,K))**2)+ &
	0.25*Ga/dr* &
	(ABS(V1(I,J,K)+V1(I,J+1,K))*(V1(I,J,K)-V1(I,J+1,K))- &
	ABS(V1(I,J-1,K)+V1(I,J,K))*(V1(I,J-1,K)-V1(I,J,K))) 
duvdth=0.25/dth*((U1(I,J,K)+U1(I,J+1,K))*(V1(I,J,K)+V1(I+1,J,K))- &
	(U1(I-1,J,K)+U1(I-1,J+1,K))*(V1(I-1,J,K)+V1(I,J,K)))+ &
	0.25*Ga/dth* &
	(ABS(U1(I,J,K)+U1(I,J+1,K))*(V1(I,J,K)-V1(I+1,J,K))- &
	ABS(U1(I-1,J,K)+U1(I-1,J+1,K))*(V1(I-1,J,K)-V1(I,J,K)))
dwvdphi=0.25/dphi*((W1(I,J,K)+W1(I,J+1,K))*(V1(I,J,K)+V1(I,J,K+1))- &
	(W1(I,J,K-1)+W1(I,J+1,K-1))*(V1(I,J,K-1)+V1(I,J,K)))+ &
	0.25*Ga/dphi* &
	(ABS(W1(I,J,K)+W1(I,J+1,K))*(V1(I,J,K)-V1(I,J,K+1))- &
	ABS(W1(I,J,K-1)+W1(I,J+1,K-1))*(V1(I,J,K-1)-V1(I,J,K)))

LV(I,J,k)= &
	D2VDr2+(2./r(j))*DVDr &
	+(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*D2VDth2 &
	+(1./(r(j)*r(j)))*D2VDphi2 &
	+(1./(r(j)*r(j)*Tan(phi(k))))*DVDphi

flr(i,j,k)=(1. - 0.5*beta*(TE2(i,j+1,k)+TE2(i,j,k)))*gr

	G(I,J,k)= &
	V1(I,J,k)+ &
	dt*( &
	(1./RE)*( &
	LV(I,J,k)-(2./(r(j)*r(j)))*V1(I,J,k)-(2./(r(j)*r(j)))*DWDphi &
	+(1./Tan(phi(k)))*W1(I,J,k) &
	-(2./(r(j)*r(j)*Sin(phi(k))))*DUDth &
	) &
	-DV2Dr-(1./(r(j)*Sin(phi(k))))*DUVDth &
	-(1./r(j))*DWVDphi &
	+(1./r(j))*(U1(I,J,k)*U1(I,J,k)+W1(I,J,k)*W1(I,J,k)-2.*V1(I,J,k)*V1(I,J,k)) &
	-(1./(r(j)*Tan(phi(k))))*(W1(I,J,k)*V1(I,J,k)) &
	+FLr(I,J,k) &
	)

enddo
enddo
enddo
return
end subroutine GCALC

subroutine LCALC(L,TE2,u1,v1,w1,flphi,r,phi,&
	imin,imax,jmin,jmax,kmin,kmax,&
	Re,Ga,dr,dth,dphi,dt,gphi,beta)
implicit none
integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: Re,Ga,dr,dth,dphi,dt,gphi,beta
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: L,TE2,u1,v1,w1,flphi
real, intent(inout) :: r(jmin-1:jmax+1),&
			phi(kmin-1:kmax+1)

integer :: i,j,k
real :: d2wdr2,dwdr,d2wdphi2,dwdphi,d2wdth2,&
	dudth,dvdphi,dvwdr,duwdth,dw2dphi

real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1) :: Lw

do i=imin,imax
do j=jmin,jmax
do k=kmin,kmax-1
d2wdr2=(w1(i,j+1,k)-2.*w1(i,j,k)+w1(i,j-1,k))/(dr*dr)
dwdr=(w1(i,j+1,k)-w1(i,j-1,k))/(2.*dr)
d2wdphi2=(w1(i,j,k+1)-2.*w1(i,j,k)+w1(i,j,k-1))/(dphi*dphi)
dwdphi=(w1(i,j,k+1) - w1(i,j,k-1))/(2.*dphi)
d2wdth2=(w1(i+1,j,k)-2.*w1(i,j,k)+w1(i-1,j,k))/(dth*dth)
dvdphi=(v1(i,j,k+1)-v1(i,j,k-1))/(2.*dphi)
dudth=(u1(i+1,j,k)-u1(i-1,j,k))/(2.*dth)
dvwdr=0.25/dr*((V1(I,J,K)+V1(I,J,K+1))*(W1(I,J,K)+W1(I,J+1,K))- &
	(V1(I,J-1,K)+V1(I,J-1,K+1))*(W1(I,J-1,K)+W1(I,J,K)))+ &
	0.25*Ga/dr* &
	(ABS(V1(I,J,K)+V1(I,J,K+1))*(W1(I,J,K)-W1(I,J+1,K))- &
	ABS(V1(I,J-1,K)+V1(I,J-1,K+1))*(W1(I,J-1,K)-W1(I,J,K)))
duwdth=0.25/dth*((U1(I,J,K)+U1(I,J,K+1))*(W1(I,J,K)+W1(I+1,J,K))- &
	(U1(I-1,J,K)+U1(I-1,J,K+1))*(W1(I-1,J,K)+W1(I,J,K)))+ &
	0.25*Ga/dth* &
	(ABS(U1(I,J,K)+U1(I,J,K+1))*(W1(I,J,K)-W1(I+1,J,K))- &
	ABS(U1(I-1,J,K)+U1(I-1,J,K+1))*(W1(I-1,J,K)-W1(I,J,K)))
dw2dphi=0.25/dphi*((W1(I,J,K)+W1(I,J,K+1))**2-(W1(I,J,K-1)+W1(I,J,K))**2)+ &
	0.25*Ga/dphi* &
	(ABS(W1(I,J,K)+W1(I,J,K+1))*(W1(I,J,K)-W1(I,J,K+1))- &
	ABS(W1(I,J,K-1)+W1(I,J,K))*(W1(I,J,K-1)-W1(I,J,K)))

LW(I,J,k)= &
	D2WDr2+(2./r(j))*DWDr &
	+(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*D2WDth2 &
	+(1./(r(j)*r(j)))*D2WDphi2 &
	+(1./(r(j)*r(j)*Tan(phi(k))))*DWDphi

flphi(i,j,k)=(1. - 0.5*beta*(TE2(i,j,k+1)+TE2(i,j,k)))*gphi

	L(I,J,k)= &
	W1(I,J,k)+ &
	dt*( &
	(1./RE)*(&
	LW(I,J,k)+(2./(r(j)*r(j)))*DVDphi &
	-(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*W1(I,J,k) &
	-(2./(r(j)*r(j)))*(Cos(phi(k))/(Sin(phi(k))*Sin(phi(k))))*DUDth &
	) &
	-DVWDr-(1./(r(j)*Sin(phi(k))))*DUWDth-(1./r(j))*DW2Dphi &
	-(3./r(j))*(V1(I,J,k)*W1(I,J,k)) &
	-(1./(r(j)*Tan(phi(k))))*(W1(I,J,k)*W1(I,J,k)-U1(I,J,k)*U1(I,J,k)) &
	+FLphi(I,J,k) &
	)
enddo
enddo
enddo
return
end subroutine LCALC

subroutine PRESSURE(P1,P2,F,G,L,r,phi,&
		imin,imax,jmin,jmax,kmin,kmax,&
		dth,dr,dphi,dt)
implicit none
integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: dth,dr,dphi,dt
real, intent(inout) :: r(jmin-1:jmax+1),phi(kmin-1:kmax+1)
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: &
	p1,p2,F,G,L
real :: HTP
integer :: i,j,k,pt
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1) :: QP

HTP=0.000001
do pt=1,100
FORALL(i=imin:imax,j=jmin:jmax,k=kmin:kmax)
Qp(i,j,k)=(1./dt)*( &
	(2./r(j))*G(I,J,k) + &
	(G(I,J,k)-G(I,J-1,k))/dr + &
	(1./(r(j)*Sin(phi(k))))*(F(I,J,k)-F(I-1,J,k))/dth + &
	(1./(r(j)*Tan(phi(k))))*L(I,J,k) + &
	(1./r(j))*(L(I,J,k)-L(I,J,k-1))/dphi &
	)

	P2(I,J,k)=P1(I,J,k)+ HTP* ( &
	(2./r(j))*(P1(I,J+1,k)-P1(I,J-1,k))/(2*dr) + &
	(P1(I,J+1,k)-2.*P1(I,J,k)+P1(I,J-1,k))/(dr**2) + &
	(1./(r(j)*r(j)*Sin(phi(k))*Sin(phi(k))))*(P1(I+1,J,k)-2.*P1(I,J,k)+P1(I-1,J,k))/(dth**2) + &
	(1./(r(j)*r(j)*Tan(phi(k))))*(P1(I,J,k+1)-P1(I,J,k-1))/(2*dphi) + &
	(1./(r(j)*r(j)))*(P1(I,J,k+1)-2.*P1(I,J,k)+P1(I,J,k-1))/(dphi**2) - &
	Qp(i,j,k) &
	)
END FORALL

FORALL(j=jmin:jmax,k=kmin:kmax)
P1(imin-1,J,k)=P1(imin,J,k)
P2(imin-1,J,k)=P2(imin,J,k)
P1(imax+1,J,k)=P1(imax,J,k)
P2(imax+1,J,k)=P2(imax,J,k)

P1(imin,J,K)=P1(imax,J,K)
P2(imin,J,k)=P2(imax,J,k)
END FORALL

FORALL(i=imin:imax,k=kmin:kmax)
P1(I,jmin-1,k)=P1(I,jmin,k)
P2(I,jmin-1,k)=P2(I,jmin,k)
P1(I,jmax+1,k)=P1(I,jmax,k)
P2(I,jmax+1,k)=P2(I,jmax,k)
END FORALL

FORALL(i=imin:imax,j=jmin:jmax)
P1(I,J,kmin-1)=P1(I,J,kmin)
P2(I,J,kmin-1)=P2(I,J,kmin)
P1(I,J,kmax+1)=P1(I,J,kmax)
P2(I,J,kmax+1)=P2(I,J,kmax)
END FORALL

P1=P2
enddo
return
end subroutine PRESSURE

subroutine U2V2W2(u2,v2,w2,F,G,L,p2,r,phi,&
		imin,imax,jmin,jmax,kmin,kmax,&
		dth,dphi,dr,dt)
implicit none
integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax
real, intent(inout) :: dth,dphi,dr,dt
real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1), intent(inout) :: &
	U2,V2,W2,F,G,L,P2
real, intent(inout) :: r(jmin-1:jmax+1),phi(kmin-1:kmax+1)
integer :: i,j,k

FORALL(i=imin:imax-1,j=jmin:jmax,k=kmin:kmax)
U2(I,J,k)=F(I,J,k)-(dt/dth)*(1./(r(j)*Sin(phi(k))))*(P2(I+1,J,k)-P2(I,J,k))
END FORALL

FORALL(i=imin:imax,j=jmin:jmax-1,k=kmin:kmax)
V2(I,J,k)=G(I,J,k)-(dt/dr)*(P2(I,J+1,k)-P2(I,J,k))
END FORALL

FORALL(i=imin:imax,j=jmin:jmax,k=kmin:kmax-1)
W2(I,J,k)=L(I,J,k)-(dt/dphi)*(1./r(j))*(P2(I,J,k+1)-P2(I,J,k))
END FORALL
return
end subroutine U2V2W2

subroutine ANIMACION(r,th,phi,u2,v2,w2,&
		imin,imax,jmin,jmax,kmin,kmax,Nplot,t,TE2)
implicit none

integer, intent(inout) :: imin,imax,jmin,jmax,kmin,kmax,Nplot,t

real, dimension(imin-1:imax+1,jmin-1:jmax+1,kmin-1:kmax+1),intent(inout) ::&
	 u2,v2,w2,TE2
real, intent(inout) :: r(jmin-1:jmax+1),phi(kmin-1:kmax+1),th(imin-1:imax+1)
integer :: i,j,k

integer :: NCOUNT
character :: EXT*4, DESTINY*512, NCTEXT*4, OUTFILE*512

EXT='.dat'
DESTINY='PATH' !Destination path, example '/garzon/Desktop/Case1/anim/', there is a limit in the characters elements
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I4.4)') NCOUNT

OUTFILE=trim(DESTINY)//'Time_'//trim(NCTEXT)//EXT

OPEN(10,file=OUTFILE)
2    format(2x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5)
write(10,*) ' TITLE = "NS2DSphere" '
write(10,*) ' VARIABLES = "X", "Y","Z", "U", "V", "W","MAGU","TEMP"'
write(10,*) ' ZONE T="',NCOUNT,'", I=',imax+2,', J=',jmax,', K=',kmax

DO k=kmin,kmax
DO j=jmin,jmax
DO i=imin-1,imax+1

!DO i=imin-1,imax+1 !For paraview

	write(10,2) r(j)*Sin(phi(k))*Cos(th(i)), &
		r(j)*Sin(phi(k))*Sin(th(i)), &
		r(j)*Cos(phi(k)), &
		V2(I,J,k)*Cos(Th(i))*Sin(phi(k))-U2(I,J,k)*Sin(Th(i))+W2(I,J,k)*Cos(Th(i))*Cos(phi(k)), &
		V2(I,J,k)*Sin(Th(i))*Sin(phi(k))+U2(I,J,k)*Cos(Th(i))+W2(I,J,k)*Sin(Th(i))*Cos(phi(k)), &
		V2(I,J,k)*Cos(phi(k))-W2(I,J,k)*Sin(phi(k)),&
		SQRT((V2(I,J,k)*Cos(Th(i))*Sin(phi(k))-U2(I,J,k)*Sin(Th(i))+W2(I,J,k)*Cos(Th(i))*Cos(phi(k)))**2&
		+(V2(I,J,k)*Sin(Th(i))*Sin(phi(k))+U2(I,J,k)*Cos(Th(i))+W2(I,J,k)*Sin(Th(i))*Cos(phi(k)))**2&
		+(V2(I,J,k)*Cos(phi(k))-W2(I,J,k)*Sin(phi(k)))**2),TE2(i,j,k)
ENDDO
ENDDO
ENDDO
return
end subroutine
