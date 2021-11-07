#!/usr/bin/env python3
import subprocess
import os
import shlex
def cl(command):
    #ip::string, command line as string input
    #op::string, return value is the output of command line
    #Notice, each time when cl is called, cl starts from current directory.
    #Use three \' if you want to input multiple lines
    #String emptiness detection is taken from https://stackoverflow.com/a/55747410
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
    (stdout_str, stderr_str) = p.communicate()
    if (not "".__eq__(stdout_str)) and (not stdout_str.isspace()):
        print(stdout_str)
    if (not "".__eq__(stderr_str)) and (not stderr_str.isspace()):
        print('The command has also written to stderr:')
        print(stderr_str)
    return stdout_str

order = input('Please input the maximum order of the polynomial: ')
symmetry = input('Please input the permutation symmetry of the molecule: ')
train_x = input('Please input the name of the data file: ')
arg = order +' '+ symmetry

print("")
print("Generating the fitting bases... (This might take time) \n")

cl('''
cd src
cd emsa
make
cp msa ../
cd ../
./msa '''+ arg +  '''
perl postemsa.pl ''' + arg + '''
perl derivative.pl ''' + arg
)

f = open(train_x)
nol = 0
for line in f:
  nol=nol+1
  if nol==1:
    natom=int(line)
nconfig = nol/(natom+2)
f.close()
print('1. Input file info:')
print(('Number of atoms is : ' + str (natom)))
print(('Number of configurations is: '+str(nconfig)+'\n'))

f = open('./src/basis.f90')
nol=1 #Num of lines in file
for line in f:
  if nol==8:
    ncoeff = int(line.split(':')[1].split(')')[0])
    ncoeff = ncoeff + 1
  if nol==24:
    nmonomial = int(line.split(':')[1].split(')')[0])
    nmonomial = nmonomial + 1
    break
  nol=nol+1
f.close()

print('2. Polynomial info:')
print(('Given polynomial order: '+ order))
print(('Given symmetry: '+ symmetry))
print(('Number of coefficients is: ' + str(ncoeff) +'\n'))

ans = input('Would you like to continue? y/n \n')
if ans == 'n' or ans == "N":
    print('Fitting program terminated')
    quit()

print("")
ans = input(
'''Do you want to include energy gradient (negative force) in the fitting? y/n
''')
if ans == 'y' or ans == 'Y':
    havegd = ".true."
else:
    havegd = ".false."

print("")
ans = input(
'''If you would like to set additional parameters, please specify the file for
these additional parameters. Otherwise enter "n" to use the default:
''')
if ans == "n" or ans == "N":
    a0 = "2.d0"
    wt = "1.e10"
else:
    f = open(ans,"r")
    line = f.readline()
    line = f.readline()
    a0 = line.split()[0]
    wt = line.split()[1]

cl('''cd src
sed 's/a = 2.0d0/a = ''' + a0 + '''/g' gradient.f90 > temp.f90
mv temp.f90 gradient.f90
''')

g = open('./src/fit.f90','w')
g.write('''program fit
use basis
use gradient
implicit none

  external dgelss
  character (len=2) :: symb
  logical :: havegrad
  real :: work(150000), dr(3), vmin, a0
  real :: rmse, grmse, wrmse, wgrmse, dwt
  integer :: i, j, k, m, info, rank
  integer :: ne, ntot, ndis, ncoeff, natm
  real,allocatable::xyz(:,:,:),v(:),g(:,:,:),wt(:)
  real,allocatable::yij(:,:),drdx(:,:,:),mono(:),poly(:),dvp(:)
  real,allocatable::A(:,:),b(:),coeff(:),s(:)
  real,allocatable::v_out(:),g_out(:,:,:)
  real,parameter::auang=0.5291772083d0
  real,parameter::aucm=219474.63d0

  havegrad=''' + havegd + '''
  natm=''' + str(natom) + '''
  ncoeff=''' + str(ncoeff)+ '''
  ne=''' + str(nconfig)+ '''
  a0=''' + a0 + '''
  dwt=''' + wt + '''

  ndis=natm*(natm-1)/2

  open(10,file="'''+train_x+'''",status="old")
  open(11,file="coeff.dat",status='unknown')
  open(12,FILE='Fit-pot.dat',status='unknown')
  open(13,file='Fit-grad.dat',status='unknown')

  allocate(xyz(ne,natm,3))
  allocate(v(ne),v_out(ne),coeff(ncoeff),s(ncoeff))
  allocate(yij(ne,ndis),wt(ne))
  allocate(poly(ncoeff))
  if (havegrad) then
     ntot=ne*(3*natm+1)
     allocate(g(ne,natm,3),g_out(ne,natm,3))
     allocate(mono('''+str(nmonomial)+'''))
     allocate(drdx(ne,3*natm,ndis))
     allocate(dvp(ncoeff))
     allocate(A(ntot,ncoeff),b(ntot))
  else
     allocate(A(ne,ncoeff),b(ne))
  end if

  do i=1,ne
     read(10,*)
     read(10,*) v(i)
     do j=1,natm
        if (havegrad) then
           read(10,*) symb,xyz(i,j,:),g(i,j,:)
        else
           read(10,*) symb,xyz(i,j,:)
        end if
     end do
  end do
  xyz = xyz / auang
  vmin = minval(v)

  yij=0.d0
  if (havegrad .eqv. .true.) drdx=0.d0
  do m=1,ne
     k = 1
     do i=1,natm-1
        do j=i+1,natm

           dr=xyz(m,i,:)-xyz(m,j,:)
           yij(m,k)=sqrt(dot_product(dr,dr))

           if (havegrad .eqv. .true.) then
              drdx(m,3*i-2:3*i,k) = dr(:)/yij(m,k)
              drdx(m,3*j-2:3*j,k) = -drdx(m,3*i-2:3*i,k)
           end if

           yij(m,k)=exp(-yij(m,k)/a0)
           k=k+1
        end do
     end do
  end do
  deallocate(xyz)

  A=0.0
  b=0.0
  do m=1,ne
     call bemsav(yij(m,:),poly)
     wt(m)=dwt/(dwt+v(m)-vmin)
     A(m,:)=poly*wt(m)
     b(m)=v(m)*wt(m)

     if (havegrad) then
        call evmono(yij(m,:),mono)
        call evpoly(mono,poly)
        do i=1,3*natm
           k=ceiling(i/3.0)
           call dbemsav(drdx(m,:,:),dvp,mono,poly,i)
           A((ne+i+(m-1)*3*natm),:)=dvp*wt(m)
           b((ne+i+(m-1)*3*natm))=g(m,k,(i-(k-1)*3))*wt(m)
        end do
     end if
  end do

  if (havegrad) then
     call dgelss(ntot,ncoeff,1,A,ntot,b,ntot,s,1.0d-11,rank,work,150000,info)
  else
     call dgelss(ne,ncoeff,1,A,ne,b,ne,s,1.0d-11,rank,work,150000,info)
  end if
  if (info .ne. 0) then
      stop 'LAPACK DGELSS solver failure'
  endif

  coeff(:)=b(1:ncoeff)
  do i=1,ncoeff
     write(11,*) coeff(i)
  end do

  write(12,"(A)") "#   V_ai (hartree)     V_PES (hartree)     Diff. (cm-1)"
  rmse=0.d0
  wrmse=0.d0
  do m=1,ne
     v_out(m)=emsav(yij(m,:),coeff)
     write (12,'(2F15.8,F12.2)') v(m),v_out(m),abs(v(m)-v_out(m))*aucm
     rmse=rmse+(v(m)-v_out(m))**2
     wrmse=wrmse+(wt(m)*(v(m)-v_out(m)))**2
  end do
  rmse=sqrt(rmse/dble(ne))
  wrmse=sqrt(wrmse/dble(ne))
  write(*,"(A,F9.2,A)") 'Overall RMSE for energy:', rmse*aucm, ' cm-1'
  write(*,"(A,F9.2,A)") "wrighted RMSE for energy:", wrmse*aucm, " cm-1"

  if (havegrad) then
     write(13,*) "#  Grad_ai (h/bohr)   Grad_pes (h/bohr)   Diff (h/bohr)"
     grmse=0.d0
     wgrmse=0.d0
     do m=1,ne
        write(13,*) natm
        call evmono(yij(m,:),mono)
        call evpoly(mono,poly)
        do i=1,3*natm
           k=ceiling(i/3.0)
           j=i-(k-1)*3
           g_out(m,k,j)=demsav(drdx(m,:,:),coeff,mono,poly,i)
           write(13,'(3F15.8)')g(m,k,j),g_out(m,k,j),abs(g(m,k,j)-g_out(m,k,j))
           grmse=grmse+(g(m,k,j)-g_out(m,k,j))**2
           wgrmse=wgrmse+(wt(m)*(g(m,k,j)-g_out(m,k,j)))**2
        end do
     end do
     grmse=sqrt(grmse/dble(3*natm*ne))
     wgrmse=sqrt(wgrmse/dble(3*natm*ne))
     write(*,"(A,F13.7,A)") 'Overall RMSE for gradient:', grmse, ' hartree/bohr'
     write(*,"(A,F13.7,A)") "Weighted RMSE for gradient:",wgrmse," hartree/bohr"
  end if

  deallocate(v,v_out,coeff,s,yij,poly,A,b)
  if (havegrad) deallocate(g,g_out,mono,dvp,drdx)

  close (10)
  close (11)
  close (12)
  close (13)
end program
''')
g.close() #Must close the file handle if you want to compile this file.
cl('''
cd src
make'''
)

print("Fitting... (This might take time) \n")

cl('''cp ./src/fit.x ./
./fit.x '''+train_x+'''
rm fit.x
mv ./src/basis.f90 ./
mv ./src/gradient.f90 ./
cp -p ./src/Makefile ./ '''
)

g = open('pes_shell.f90','w')
g.write('''module pes_shell
use basis
use gradient
implicit none

  real::coeff(1:'''+str(ncoeff)+''') ! change to number of coefficients
                     ! (size of c in bemsa.f90)
  save coeff

contains
  !==================================!
  ! read the coefficients of the PES !
  !==================================!
  subroutine pes_init()
    !::::::::::::::::::
    integer::i

    open(10,file='coeff.dat',status='old')

    do i=1,size(coeff)
       read (10,*) coeff(i)
    end do

    return
    close (10)
  end subroutine pes_init

  !====================================!
  ! Function to evaluate the potential !
  !====================================!
  function f(xyz)
    real,dimension(:,:),intent(in)::xyz
    real::f
    !::::::::::::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::x
    real,dimension(3)::dr
    real::a0  ! the same as the fitting code
    integer::i,j,k

    a0 = ''' + a0 + '''

    k = 1
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr = xyz(:,i) - xyz(:,j)
          x(k) = sqrt(dot_product(dr,dr))
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-x(i)/a0)
    end do

    f=emsav(x,coeff)

    return
  end function f

  !===========================!
  ! function to calculate the !
  ! analytical gradient       !
  !===========================!
  function g(xyz)
    real,dimension(:,:),intent(in)::xyz
    real,dimension(size(xyz,2)*3)::g
    !::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::r,x
    real,dimension(3,size(xyz,2)*(size(xyz,2)-1)/2)::dr
    real,dimension(size(xyz,2)*3,size(xyz,2)*(size(xyz,2)-1)/2)::drdx
    real,dimension(1:''' + str(ncoeff) + ''')::p   ! change to number of popynomials
                               ! (size of p in bemsa.f90)
    real,dimension(1:'''+str(nmonomial)+''')::m  ! change to number of monomials
                              ! (size of m in bemsa.f90
    real::a0
    integer::i,j,k

    a0 = ''' + a0 + '''

    k = 1
    drdx = 0.d0
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr(:,k) = xyz(:,i) - xyz(:,j)
          r(k) = sqrt(dot_product(dr(:,k),dr(:,k)))

          drdx(3*i-2:3*i,k) = dr(:,k)/r(k)
          drdx(3*j-2:3*j,k) = -drdx(3*i-2:3*i,k)
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-r(i)/a0)
    end do

    call evmono(x,m)
    call evpoly(m,p)

    do i=1,3*size(xyz,2)
       g(i) = demsav(drdx,coeff,m,p,i)
    end do

    return
  end function g

end module pes_shell
''')
g.close()

g = open("test.f90","w")
g.write('''
program main
use pes_shell
implicit none

  real :: v, v0
  logical :: havegrad
  character (len=2) :: symb
  integer :: i, natm
  real, dimension (:,:), allocatable :: xyz
  real, dimension (:), allocatable :: grad, g0
  character(len=32)::fname, bname
  real,parameter::auang= 0.5291772083d0

  havegrad = '''+havegd+'''

  call getarg(1,fname)
  i=index(fname,'.',.true.)
  bname=fname(1:i-1)

  open(12,file=trim(fname),status='old')
  open(13,file=trim(bname)//'.out',status='unknown')

  call pes_init()
  read (12,*) natm
  read (12,*) v0

  allocate(xyz(3,natm))
  if (havegrad) allocate(grad(3*natm),g0(3*natm))

  do i=1,natm
     if (havegrad) then
        read (12,*)symb, xyz(:,i), g0(3*i-2:3*i)
     else
        read (12,*)symb, xyz(:,i)
     end if
  end do

  xyz=xyz/auang
  v=f(xyz)
  write(13,'(A)') "# E_ai (hartree)   E_fit (hartree)"
  write(13,'(2F16.8)') v0, v

  if (havegrad) then
     grad=g(xyz)
     write(13,'(A)') "# Grad_ai (hartree/bohr)  Grad_fit:"
     do i=1,size(grad)
        write(13,'(2F16.8)') g0(i), grad(i)
     end do
     deallocate(g0,grad)
  end if

  deallocate(xyz)
  close (12)
  close (13)

end program
''')
g.close()

cl('''make test.x
cp ./src/test.xyz ./'''
)

print ('4. In order to run the test program, use command:')
print ('./test.x test.xyz \n')
print ('End of program')
