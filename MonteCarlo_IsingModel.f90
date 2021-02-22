!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                        !!!!!!
!!!!! Cristina Outeda Rúa                                         22/12/2020 !!!!!!
!!!!!                        Exercise 2D Ising Model                         !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program modelo_ising_2d
    implicit none
    integer :: N , L, D, E , M , k , Z , iterations, counter
    real(4), dimension (:), allocatable :: table
    integer, dimension (:), allocatable :: s , E_sist
    integer, dimension (:,:), allocatable :: nbr
    double precision ::  T , kb
    real(4) :: beta 
    real(4) :: r1279
    external :: toroidal_boundary_cond , neighbours, observables
    call setr1279(3497298)
    D = 2
    Z = 2 * D
    L = 4
    N = L ** 2
    T = 2.0d0
    kb = 1.0d0   
    beta = real(1.d0 / (kb * T))
    !print*, beta
    iterations = 10 ** 6
    allocate(s(N),nbr(Z,N),table(-2*Z:2*Z),E_sist(iterations))
    call tabulate(Z,beta,table)
    do k = 1 , N
        s(k) = 2 * mod(nint(2.0*r1279()),2) -1
    end do
    counter = 0
    open(unit = 11, file = "ising2d.txt", action = "write")
    do k = 1 , iterations
        print*,real(k)*100/real(iterations),"%"
        call rnd_spin(L,Z,s,table)
        call observables(L,Z,s,E,M)
        write(unit=11,fmt=*)E
        counter = counter + 1
        !print*, counter, "out of" , iterations, "Percent", (real(counter)/real(iterations))*100.
        E_sist(k) = E
    end do
    close(11)
    call binning(E_sist,iterations)
end program modelo_ising_2d

subroutine binning(E_1, iterations)
    integer :: iterations , j , k, counter, bin
    integer :: E_1(iterations)
    real(4), dimension (:,:), allocatable :: E_med
    real(4), dimension (:), allocatable :: E_var
    real(4) :: average
    bin = 10
    allocate(E_med(iterations,bin))
    allocate(E_var(bin))
    E_med = 0.d0
    E_var = 0.d0
    E = 0.d0
    average = real(sum(E_1(:)))/real(size(E_1))
    !print*, average
    do j = 1 , bin
    do k = 2, iterations
        if (mod(k,2**j)==0) then
            E_med(k/2**j,j) = real(sum(E_1(k-2**j+1:k))) / real(2**j)
        end if
    end do
    counter = 0     
    do k = 1, iterations
        if (abs(E_med(k,j))>10d0**(-6)) then
            E_var(j) =  E_var(j) + ((average - E_med(k,j))**2)
            counter = counter + 1
        end if
        end do
    E_var(j) = E_var(j) / (real(counter-1) * real(counter))
    !print*, E_var(j) , (real(counter-1) * real(counter))
end do
open(unit=12,file="variance.dat",action="write")
do k = 1, bin
write(12,fmt=*)2**k , E_var(k)
end do
close(12)
end subroutine


subroutine rnd_spin(L,Z,s,table)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!! Metropolis without proper random
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer :: i , Z, L , dE, N , k
integer :: S(L**2), saux(L**2) , h(L**2) , nbr(Z,L**2)
real(4) :: table(-2*Z:2*Z) , r1279
call setr1279(43982)
N = L ** 2
do k = 1 , N
    i = mod(nint(real(N)*r1279()) ,N) + 1
    saux(i) = 2*mod(nint(2.0*r1279()),2) -1
    !print*, i, saux(i)
    call neighbours(Z,L,nbr)                                           
    h(i) = sum(s(nbr(:,i)))
    dE = 2 * saux(i) * h(i)
    if (dE.lt.0.or.dE.gt.0.and.r1279().lt.table(dE)) s(i) = saux(i)
end do
end subroutine

subroutine tabulate(Z,beta,table)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!! Table of probabilities for dE
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer :: Z , dE
    real(4) ::  beta, table(-2*Z:2*Z)
    do dE = -2*Z , 2*Z , 2
        table(dE) = exp(-beta* real(dE))
    end do
end subroutine

subroutine toroidal_boundary_cond(L,in)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!
    !!!!! Periodic Toroidal Boundary Condidions
    !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none                                                        !!!!!!!
    integer :: i, L , in(2,L)                                            !!!!!!!
    do i = 1, L                                                          !!!!!!!
        in(1,i)=i-1                                                      !!!!!!!
        in(2,i)=i+1                                                      !!!!!!!
        enddo                                                            !!!!!!!  
    in(1,1)=L                                                            !!!!!!!
    in(2,L)=1                                                            !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

subroutine neighbours(Z,L,nbr)
    integer :: Z , L, in(2,L) , nbr(Z,L**2)
    integer :: i , x , y 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! computation of the neighbours for a square lattice
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call toroidal_boundary_cond(L,in)
i=0                                                                      !!!!!!!
do y = 1, L                                                              !!!!!!!
    do x = 1, L                                                          !!!!!!!
        i=i+1                                                            !!!!!!!
        nbr(1,i) = in(2,x) + L*(y-1) ! right                             !!!!!!!
        nbr(2,i) = in(1,x) + L*(y-1) ! left                              !!!!!!!
        nbr(3,i) = x + L*(in(2,y)-1) ! up                                !!!!!!!
        nbr(4,i) = x + L*(in(1,y)-1) ! down                              !!!!!!!
    enddo                                                               !!!!!!!
enddo                                                                   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine                                                                                              

subroutine observables(L,Z,s,E,M)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!  Read the input file and fill the vector with the spin values
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none                                                       !!!!!!!!!!
    integer :: M , E, Eint, Z, L , nbr(Z,L**2), s(L**2)                 !!!!!!!!!!
    integer :: i , k                                                    !!!!!!!!!!
    call neighbours(Z,L,nbr)
    E=0                                                                 !!!!!!!!!!
    M=0                                                                 !!!!!!!!!!
    do i=1,L**2                                                         !!!!!!!!!!
        Eint=0                                                          !!!!!!!!!!
        do k=1,Z                                                        !!!!!!!!!!
            Eint=Eint+s(nbr(k,i))                                       !!!!!!!!!!
        enddo                                                           !!!!!!!!!!
        E = E + Eint*s(i)/2                                             !!!!!!!!!!
        M = M + s(i)                                                    !!!!!!!!!!
    enddo                                                               !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    print*,"Energy: ",E,"Magnetitzation: ", M                          !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine                                                       
