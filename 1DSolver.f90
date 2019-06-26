!this program uses Guassian-Lobatto quadrature and the DVR method to solve the schrodinger equation in 3D. 
!Uses sparce matricies for T and V
!Outputs can be eigenvalues  
program Solver
!use mkl_spblas
implicit none
        real*8, allocatable :: nodesX(:), weightsX(:), nodesY(:), weightsY(:), nodesZ(:), &
                weightsZ(:),holder(:),EigenVals(:), res(:), Hsparse(:)
        real*8, allocatable :: Xx(:,:),dXx(:,:),Xy(:,:),dXy(:,:),Xz(:,:),dXz(:,:), &
                weightsDMX(:,:),weightsDMY(:,:),weightsDMZ(:,:),EigenVecs(:,:), &
                TmatX(:,:),TmatY(:,:),TmatZ(:,:),dXxtrim(:,:), dXytrim(:,:), dXztrim(:,:), &
                dfGLX(:,:),dfGLY(:,:),dfGLZ(:,:),dfGLXt(:,:),dfGLYt(:,:),dfGLZt(:,:)
        real*8, allocatable :: indexOf(:,:,:)
        real*8 :: n,w,valX,valXp,a,b,epsout,temp,mup,mdn,r0,d,mu,dm
        integer :: feastparam(64),x,i,j,k,ip,jp,kp,m,sX,sY,sZ,sMax,numN0,ijk,ijkp,sXc,sYc,sZc, &
                info,loop,rCount,row,m0,fsx,fsy,fsz,fsmax,fxmax
        integer, allocatable :: Hrow(:), Hcol(:)
        integer :: Tstart, Tend, rate
        real*8, external :: x12,x13,x14,x23,x24,x34
 
        open(unit=1,file="nodesAndWeights10.dat")
        open(unit=2,file="nodesAndWeights20.dat")
        open(unit=3,file="nodesAndWeights30.dat")
        open(unit=4,file="nodesAndWeights40.dat")
        open(unit=5,file="nodesAndWeights50.dat")
        open(unit=6,file="nodesAndWeights60.dat")
        !open(unit=7,file="nodesAndWeights70.dat")
        !open(unit=8,file="nodesAndWeights80.dat")
        !open(unit=9,file="nodesAndWeights90.dat")
        !open(unit=10,file="nodesAndWeights100.dat")
        
        open(unit=100,file="dataOut.dat")

        !write(*,*) "enter rescale range"
        !read(*,*) a,b 
        !a=-5d0
        !b=5d0
        !do x=1,6
        x=1

        print *, "starting set",x
        call system_clock(Tstart)

        read(x,*) sX
        sY=sX
        sZ=sX
        sXc=sX-2
        sYc=sY-2
        sZc=sZ-2
        sMax=sXc*sYc*sZc


        fsx=sx/2-1
        fsy=sy/2-1
        fsz=sz/2-1
        fsMax=fsx*fsy*fsz

        fxMax=sXc*sYc*fsz
        
        dm=5d0
        mup=1d0
        mdn=1d0
        r0=0.2d0
        d=0d0
        a=1/r0
        !mu=(1d0/4d0)**(-1d0/3d0)
        !mu = 1d0

        allocate(nodesX(sX),weightsX(sX),nodesY(sY),weightsY(sY),nodesZ(sZ),weightsZ(sZ),holder(sX), &
                Xx(sX,sX),dXx(sX,sX),Xy(sY,sY),dXy(sY,sY),Xz(sZ,sZ),dXz(sZ,sZ), &
                weightsDMX(sX,sX),weightsDMY(sY,sY),weightsDMZ(sZ/2,sZ/2),indexOf(sXc,sYc,fsz), &
                TmatX(sXc,sXc), TmatY(sYc,sYc),TmatZ(fsz,fsz),dXxtrim(sXc,sX),dXytrim(sYc,sY),dXztrim(sZc,sZ), &
                dfGLX(fsx,sx),dFGLY(fsy,sy),dfGLZ(fsz,sz),dfGLXt(fsx,sx/2),dfGLYt(fsy,sy/2),dfGLZt(fsz,sz/2))

        !read in data for each basis
        do i=1,sX
                read(x,*) n,w
                nodesX(i)=n
                weightsX(i)=w
        end do
        nodesY=nodesX
        weightsY=weightsX
        nodesZ=nodesX
        weightsZ=weightsX
        
        !rescale the basis functions to the desired range
        call rescale(sX,0d0,dm,nodesX)
        call rescale(sY,0d0,dm,nodesY)
        call rescale(sZ,-dm,dm,nodesZ)
        !write(*,*) "rescaled"
        
        !build the DVR basis functions and their derivatives from the Gausian Lobatto nodes for each dimention
        do i=1,sX
                do j=1,sX
                        call DVRFunctions(holder,nodesX,nodesX(j),i,sX,valX,weightsX)
                        Xx(i,j)=valX

                        call DifferentiateChi(holder,nodesX,nodesX(j),i,sX,valXp,weightsX)
                        dXx(i,j)=valXp
                end do
        end do
        do i=1,sY
                do j=1,sY
                        call DVRFunctions(holder,nodesY,nodesY(j),i,sY,valX,weightsY)
                        Xy(i,j)=valX

                        call DifferentiateChi(holder,nodesY,nodesY(j),i,sY,valXp,weightsY)
                        dXy(i,j)=valXp
                end do
        end do
        do i=1,sZ
                do j=1,sZ
                        call DVRFunctions(holder,nodesZ,nodesZ(j),i,sZ,valX,weightsZ)
                        Xz(i,j)=valX

                        call DifferentiateChi(holder,nodesZ,nodesZ(j),i,sZ,valXp,weightsZ)
                        dXz(i,j)=valXp
                end do
        end do

        !print *, "got Chi's"

        !derivatives must be resized due to end points being non-continuous (imposed to zero)
        dXxtrim = dXx(2:sX-1,1:sX)
        dXytrim = dXy(2:sY-1,1:sY)
        dXztrim = dXz(2:sZ-1,1:sZ)

        print *, size(dxxtrim)

        !call xTOf(sXc,fsx,XGLXtrim,fGLX)
        !call xTOf(sXc,fsx,dXXtrim,dfGLX)
        !call xTOf(sYc,fsy,dXYtrim,dfGLY)
        call xTOf(sZc,fsz,dXZtrim,dfGLZ)
        
        !dfGLXt = dfGLX(1:fsx,1:sX/2)
        !dfGLYt = dfGLY(1:fsy,1:sY/2)
        dfGLZt = dfGLZ(1:fsz,1:sZ/2)

        !make index array to help build Tsparse matrix
        row = 0
        do i=1,sXc
                do j=1,sYc
                        do k=1,fsz
                                row=row+1
                                indexOf(i,j,k)=row
                        end do
                end do
        end do

        !make T kenetic evergy matrix for each dimenstion
        !first make diagonal weight matrices
        weightsDMX=0d0
        weightsDMY=0d0
        weightsDMZ=0d0
        do i=1,sX
                do ip=1,sX
                        if (i == ip) then
                                weightsDMX(i,ip)=weightsX(i)
                        end if
                end do
        end do
        do i=1,sY
                do ip=1,sY
                        if (i == ip) then
                                weightsDMY(i,ip)=weightsY(i)
                        end if
                end do
        end do
        do i=1,sZ/2
                do ip=1,sZ/2
                        if (i == ip) then
                                weightsDMZ(i,ip)=2d0*weightsZ(i)
                        end if
                end do
        end do
        !now make the Tmat for each dimention 
        !- this formula comes from integration by parts, and the surface term goes to zero
        TmatX = (0.5d0)*MATMUL(dXxtrim,MATMUL(weightsDMX,TRANSPOSE(dXXtrim)))
        TmatY = (0.5d0)*MATMUL(dXytrim,MATMUL(weightsDMY,TRANSPOSE(dXYtrim)))
        TmatZ = (0.5d0)*MATMUL(dfGLZt,MATMUL(weightsDMZ,TRANSPOSE(dfGLZt)))
        
        !now combine them using delta properties, and construct the Hsparse matrices explicitly, 
        !skipping creating T or V matrices
        !first loop is a fake to find array sizes to allocate memory, then second loop makes H
        numN0=0
        row=0
        rCount=1
        do m=1,2
                do i=1,sXc
                        do j=1,sYc
                                do k=1,fsz
                                        row=row+1
                                        do ip=1,sXc
                                                do jp=1,sYc
                                                        do kp=1,fsz
                                                                ijk=indexOf(i,j,k)
                                                                ijkp=indexOf(ip,jp,kp)
                                                                temp = 0d0
                                                                if((j==jp).and.(k==kp)) then
                                                                        temp=temp+TmatX(i,ip)
                                                                end if
                                                                if((i==ip).and.(k==kp)) then
                                                                        temp=temp+TmatY(j,jp)
                                                                end if
                                                                if((j==jp).and.(i==ip)) then
                                                                        temp=temp+TmatZ(k,kp)
                                                                end if
                                                                if((i==ip).and.(j==jp).and.(k==kp)) then
                                                                        Call V3d(nodesX(i+1),nodesY(j+1),nodesZ(k+1),valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x12(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x13(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x14(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x23(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x24(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                        call V2bMorse(x34(nodesX(i+1),nodesY(j+1),nodesZ(k+1),mup,mdn),r0,d,a,valX)
                                                                        temp=temp+valX
                                                                end if
                                                                if(temp /= 0) then
                                                                        numN0=numN0+1
                                                                        rCount=rCount+1
                                                                        if(m==2) then
                                                                                Hsparse(numN0)=temp
                                                                                Hcol(numN0)=ijkp
                                                                                Hrow(ijk+1)=Hrow(ijk+1)+1
                                                                        end if
                                                                end if
                                                        end do
                                                end do
                                        end do
                                end do
                        end do
                end do
                if(m==1) then
                        allocate(Hsparse(1:numN0),Hcol(1:numN0),Hrow(1:fxMax+2))
                        Hrow=0
                        Hrow(1)=1
                        numN0=0
                        row=0
                        rCount=1
                end if
        end do
        !sum up row counts for Hrow
        do i=2,fxMax+1
                Hrow(i)=Hrow(i)+Hrow(i-1)
        end do

        !print *, "Got H"
        print *, "If",size(Hsparse)," =",Hrow(fxMax+1)-1," =",size(Hcol)," then good"

        
        !set up solver and call solver
        if(x==1) then
                m0=100
        end if
        if(x==2) then
                m0=200
        end if
        if(x==3) then
                m0=300
        end if
        if(x==4) then
                m0=400
        end if
        if(x==5) then
                m0=500
        end if
        if(x==6) then
                m0=600
        end if
        allocate(EigenVals(m0),EigenVecs(fxMax,m0),res(m0))
        call feastinit(feastparam)
        feastparam(1)=1
        feastparam(2)=20
        !feastparam(4)=3
        feastparam(17)=0
        call dfeast_scsrev('F',fxMax,Hsparse,Hrow,Hcol,feastparam,epsout,loop,-6d0*d,10d0,m0,EigenVals,EigenVecs,m,res,info)
        
        print *, "info: ",info
        call system_clock(Tend,rate)
        write(100,*) x,(Tend-Tstart)/rate
        write(100,*) "" 
        
        write(100,*) (EigenVals(1)-3.5d0)/3.5d0
        write(100,*) (EigenVals(2)-5.5d0)/5.5d0
        write(100,*) (EigenVals(5)-7.5d0)/7.5d0
        write(100,*) (EigenVals(11)-9.5d0)/8.5d0

        do i=1,20
                print *, eigenvals(i)
        end do
        deallocate(nodesX,weightsX,nodesY,weightsY,nodesZ,weightsZ,holder, &
                Xx,dXx,Xy,dXy,Xz,dXz, &
                weightsDMX,weightsDMY,weightsDMZ,indexOf, &
                TmatX, TmatY,TmatZ,dXxtrim,dXytrim,dXztrim, &
                dfGLX,dFGLY,dfGLZ,dfGLXt,dfGLYt,dfGLZt)
        deallocate(Hsparse,Hcol,Hrow)
        deallocate(EigenVals,EigenVecs,res)
        
        !end do

        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        !close(7)
        !close(8)
        !close(9)
        !close(10)
        close(100)
end program Solver


!rescales the nodes of size n from a to b
subroutine rescale(n,a,b,nodes)
implicit none
        integer :: n,i
        real*8 :: nodes(*)
        real*8 :: x1,xN,a,b
        x1= nodes(1)
        xN= nodes(n)
        do i=1,n
                if((xN-x1)*(b-a)==0) then
                        nodes(i)=0d0
                else 
                        nodes(i)=a+(nodes(i)-x1)/(xN-x1)*(b-a)
                end if
        end do
end subroutine rescale


!Does the lobatto product into array
subroutine DVRFunctions(holder,nodes,node,i,sz,val,weights)
implicit none
        integer :: i,j,sz
        real*8 :: nodes(sz),holder(sz),weights(sz)
        real*8 :: val,node

        do j=1,sz
                if((nodes(i)-nodes(j))==0d0) then
                    holder(j)=0d0
                else 
                    holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
                end if
        end do
        holder(i)=1d0
        val=PRODUCT(holder)*(weights(i)**(-.5d0))
end subroutine DVRFunctions
                

!find derivitave of the nodes
subroutine DifferentiateChi(holder,nodes,node,i,sz,Xsum,weights)
implicit none
        integer :: i,sz,j,k
        real*8 nodes(sz),holder(sz),weights(sz)
        real*8 val,node,Xsum

        Xsum=0
        do k=1,sz
                if(k /= i) then
                        do j=1,sz
                                if((nodes(i)-nodes(j))==0d0) then
                                        holder(j)=0d0
                                else 
                                        holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
                                end if
                        end do
                        holder(i)=1d0
                        holder(k)=1d0
                        if((nodes(i)-nodes(k))==0) then
                                val=0d0
                        else 
                                val=PRODUCT(holder)/((nodes(i)-nodes(k)))
                        end if
                        Xsum = Xsum + val
                end if
        end do
        Xsum=Xsum*(weights(i)**(-.5d0))
end subroutine DifferentiateChi

subroutine xTOf(sz,fs,chi,f)
implicit none
        integer :: i,j,sz,fs
        real*8 chi(sz,sz+2),f(fs,sz+2)
        do i=1,fs
               do j=1,sz+2
                        f(i,j)=(1d0/sqrt(2d0))*(chi(i,j)+chi(sz+1-i,j))
                end do
        end do
end subroutine

subroutine VPotHO(x,res)
implicit none
        real*8 x,res
        
        res = (0.5d0)*x**2
end subroutine VPotHO

subroutine Vsech2(x,res)
implicit none
        real*8 x,res

        res = -((1d0)/cosh(x))**2
end subroutine Vsech2

subroutine V3D(x,y,z,res)
implicit none
        real*8 x,y,z,res
        
        res = (0.5d0)*(x**2+y**2+z**2)
end subroutine V3D


function x12(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res=(2**0.3333333333333333*mdn**0.8333333333333334*mup*x)/&
                (Sqrt((mdn*mup)/(mdn + mup))*(mup*(mdn + mup))**0.6666666666666666)
end function x12

function x13(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res = (mup**1.6666666666666667*(mdn**7*(mdn + mup))**0.16666666666666666*x + &
                mup**0.6666666666666666*(mdn**13*(mdn + mup))**0.16666666666666666*x - &
                mdn**1.6666666666666667*(mup**7*(mdn + mup))**0.16666666666666666*y - &
                mdn**0.6666666666666666*(mup**13*(mdn + mup))**0.16666666666666666*y + &
                mdn*(mdn*mup*(mdn + mup))**0.6666666666666666*z + mup*(mdn*mup*(mdn + &
                mup))**0.6666666666666666*z)/(2**0.6666666666666666*(mdn*mup)**0.8333333333333334*(mdn + &
                mup)**1.3333333333333333)
end function x13

function x14(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res = (mup**1.6666666666666667*(mdn**7*(mdn + mup))**0.16666666666666666*x + &
                mup**0.6666666666666666*(mdn**13*(mdn + mup))**0.16666666666666666*x + &
                mdn**1.6666666666666667*(mup**7*(mdn + mup))**0.16666666666666666*y + &
                mdn**0.6666666666666666*(mup**13*(mdn + mup))**0.16666666666666666*y + &
                mdn*(mdn*mup*(mdn + mup))**0.6666666666666666*z + &
                mup*(mdn*mup*(mdn + mup))**0.6666666666666666*z)/&
                (2**0.6666666666666666*(mdn*mup)**0.8333333333333334*(mdn + mup)**1.3333333333333333)
end function x14

function x23(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res = -((mup**1.6666666666666667*(mdn**7*(mdn + mup))**0.16666666666666666*x + &
                mup**0.6666666666666666*(mdn**13*(mdn + mup))**0.16666666666666666*x + &
                mdn**1.6666666666666667*(mup**7*(mdn + mup))**0.16666666666666666*y + &
                mdn**0.6666666666666666*(mup**13*(mdn + mup))**0.16666666666666666*y - &
                mdn*(mdn*mup*(mdn + mup))**0.6666666666666666*z - &
                mup*(mdn*mup*(mdn + mup))**0.6666666666666666*z)/&
                (2**0.6666666666666666*(mdn*mup)**0.8333333333333334*(mdn + mup)**1.3333333333333333))
end function x23

function x24(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res = (-(mup**1.6666666666666667*(mdn**7*(mdn + mup))**0.16666666666666666*x) - &
                mup**0.6666666666666666*(mdn**13*(mdn + mup))**0.16666666666666666*x + &
                mdn**1.6666666666666667*(mup**7*(mdn + mup))**0.16666666666666666*y + &
                mdn**0.6666666666666666*(mup**13*(mdn + mup))**0.16666666666666666*y + &
                mdn*(mdn*mup*(mdn + mup))**0.6666666666666666*z + &
                mup*(mdn*mup*(mdn + mup))**0.6666666666666666*z)/&
                (2**0.6666666666666666*(mdn*mup)**0.8333333333333334*(mdn + mup)**1.3333333333333333)
end function x24
           
function x34(x,y,z,mup,mdn) result(res)
implicit none
        real*8 x,y,z,mup,mdn,res
        res = (2**0.3333333333333333*mup**0.3333333333333333*y)/(mdn*(mdn + mup))**0.16666666666666666
end function x34

subroutine v2bMorse(r,r0,d,a,res)
implicit none
        real*8 x,res,r,r0,d,a
        res=d*(1d0-exp(-a*(abs(r)-r0)))**2d0-d
end subroutine v2bMorse

