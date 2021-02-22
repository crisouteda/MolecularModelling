program final
    implicit none
    double precision :: ro , lbox , dt, kin , mom , pot , T_instant , preasure, density(5)
    integer :: particles, i , geom , k , p
    double precision, dimension (:,:), allocatable :: r,v,fu
    external ::  boundary , lennardjones , andersen_termo , kinetic
    character(len=100)  :: filename
    particles = 130 ! N > 100
    geom = nint(dble(particles)**(1.d0/3.d0))
    particles = geom ** 3
    ro = 0.8d0
    allocate(r(particles,3))
    allocate(v(particles,3))
    allocate(fu(particles,3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                   Exercise A
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt = 0.0001 !set in 10 ^ (-4)
    call system_sc(ro,r,particles,lbox)
    call write_configuration(r,particles,"bfore_mel.xyz")
    do i = 1 , 10000
        call verlet(particles,dt,r,v,fu,lbox, pot, preasure) 
        call andersen_termo(v,particles,1000.d0)
        call kinetic(v,particles,kin,mom)
        call temperature(v,particles,T_instant)
    end do
    call write_configuration(r,particles,"after_mel.xyz")
    call inizalize_velocities(v,particles,1.5d0)
    open (unit=11,file="initial_thermodynamics.xyz",action="write")
    write(11,fmt=*)"t        , T_instant        ,kin         ,mom         ,pot       "
    dt = 0.0001
    do i = 1 , 50000 
        call verlet(particles,dt,r,v,fu,lbox, pot, preasure) 
        call kinetic(v,particles,kin,mom)
        call temperature(v,particles,T_instant)
        write(11,fmt=*)dble(i)*dt , T_instant , kin , mom , pot, preasure
    end do
    call write_configuration(r,particles,"after_iso.xyz")
    close(11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                   Exercise B
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dt = 0.0001d0
density = (/0.1d0 , 0.2d0, 0.4d0, 0.6d0, 0.8d0/)
do k = 1, 5
        ro = density(k)
        call system_sc(ro,r,particles,lbox)
        do i = 1 , 10000
            call verlet(particles,dt,r,v,fu,lbox, pot, preasure) 
            call andersen_termo(v,particles,1000.d0)
            call kinetic(v,particles,kin,mom)
            call temperature(v,particles,T_instant)
        end do

        write(filename , *) "configuration",real(density(k))
        open (unit=10,file=filename,action="write")
            write(10,fmt=*)particles
            write(10,fmt=*) 
                do p = 1 , particles
                    write(10,fmt=*)"Xe",r(p,1),r(p,2),r(p,3)
                end do
        close(10)

        write(filename , *)"termodynamics",real(density(k))
        open (unit=30,file=filename,action="write")
        !write(30,fmt=*) particles, "density = ", ro
        !write(30,fmt=*)
            do i = 1 , 50000
                call verlet(particles,dt,r,v,fu,lbox,pot, preasure)
                call andersen_termo(v,particles,1.5d0)
                call kinetic(v,particles,kin,mom)
                call temperature(v,particles,T_instant)
                write(30,fmt=*)dble(i)*dt,",",T_instant,",",kin,",",mom,",",pot,",",preasure+T_instant*ro
            end do
        close (unit=30)
 end do

end program final



subroutine system_sc(ro,r,particles,lbox)
    implicit none
    integer :: M , x, y, z, p , particles
    double precision :: ro , lbox , a, r(particles,3)
    M = nint(dble(particles)**(1.d0/3.d0))
    particles = M ** 3
    a = (1.d0/ro)**(1.0d0/3.0d0)
    lbox = a * dble(M)
    p = 0
    do x = 0, M -1
        do y = 0, M -1 
            do z = 0, M -1
                p = p + 1
                r(p,1) = dble(x)*a 
                r(p,2) = dble(y)*a
                r(p,3) = dble(z)*a
            end do
        end do
    end do
end subroutine

subroutine inizalize_velocities(v,particles,T)
    implicit none
    integer :: particles , p , c
    double precision :: v(particles,3) , T , kin
    kin = 0.d0
    do p = 1 , particles
        do c = 1, 3
            v(p,c) = rand() - 0.5d0
            kin = kin + (1.d0/2.d0) * v(p,c) ** 2.d0
        end do
    end do
    do p = 1 , particles
        do c = 1, 3
            v(p,c) = v(p,c) * sqrt(dble(particles)*T/(2.d0*kin)) 
        end do
    end do
  end subroutine

subroutine andersen_termo(v,particles,T) ! andersen thermostat
    implicit none
    integer :: particles , p , c
    double precision :: v(particles,3) , nu
    double precision :: T, PI , sigma, X1, X2, X3
    nu = 0.1d0
    PI = 4.d0*datan(1.d0)
    sigma = dsqrt(T)
    do p= 1,particles
        x1 = rand()
        if (x1<nu) then
            do c = 1, 3
                x2 = rand()
                x3 = rand()
                v(p,c)= sigma*dsqrt(-2.d0*(dlog(1.d0-x2)))*dcos(2.d0*PI*x3)
            end do
        end if
    end do
  end subroutine


  subroutine lennardjones(particles,r,fu,lbox,pot,preasure)
    implicit none
    integer :: p, c, q, particles, counter
    double precision, dimension (:,:) :: r(particles,3), fu(particles,3)
    double precision :: d(3), lbox, pot, d2 ,d4 , d6 ,d8 , d12, d14, cf6, cf12, preasure
    d = 0.d0
    fu = 0.d0
    pot = 0.d0
    preasure = 0.d0
    counter = 0
    do p = 1 , particles
    do q = 1 , particles
        if (p.ne.q) then
            do c = 1 , 3
                d(c) = r(p,c)-r(q,c)
                call boundary(d(c),lbox)
            end do
                d2 = (d(1)**2 + d(2)**2 + d(3)**2)
            if (d2.lt.(lbox/2.d0)**2.and.d2.ne.0.d0) then
                counter=counter+1
                d4 = d2 * d2
                d6 = d2 * d4
                d8 = d2 * d6
                d12 = d6 * d6
                d14 = d12 * d2
                cf6 = (lbox/2.d0)**6d0
                cf12 = cf6*cf6
                fu(p,:) = fu(p,:) + (48.d0/d14-24.d0/d8)*d(:)
                preasure = preasure + 2.d0*(1.d0/D12 - 1.d0/D6) - 2.d0*(1.d0/cf12 - 1.d0/cf6)*sqrt(d2)
                pot = pot + 2.d0*(1.d0/D12 - 1.d0/D6) - 2.d0*(1.d0/cf12 - 1.d0/cf6)  
            end if
        end if
end do
end do
preasure=((preasure/dble(counter))/(3.d0*lbox**3.d0))
end subroutine


subroutine verlet(particles,dt,r,v,fu,lbox,e_pot, preasure) 
    implicit none
    integer :: particles, p, c
    double precision :: dt ,lbox, pot , e_pot, preasure
    double precision, dimension (:,:) :: r(particles,3), v(particles,3), fu(particles,3)
    e_pot = 0.d0
    do p = 1, particles
        do c = 1 , 3
            r(p,c) = r(p,c) + v(p,c) * dt + fu(p,c) * dt * dt * 0.5d0
            call boundary(r(p,c),lbox)
            v(p,c) = v(p,c) + fu(p,c)* dt / 2.d0
        end do
    end do
    call lennardjones(particles,r,fu,lbox,pot,preasure)
    e_pot = e_pot + pot
    do p = 1, particles
        do c = 1, 3
            v(p,c) = v(p,c) + fu(p,c)* dt / 2.d0
    end do
    end do
end subroutine

subroutine boundary(x,lbox)
    implicit none
    double precision :: lbox , x
     if (x > lbox/2.0d0) x = x - lbox
     if (x < (-lbox/2.0d0)) x = x + lbox
end subroutine

subroutine write_configuration(r,particles,filename)
    implicit none
    integer :: p , particles
    double precision :: r(particles,3)
    character (len = 13) :: filename
    open(unit=10,file=filename, action = "write")
    write(10,fmt=*)particles
    write(10,fmt=*) 
        do p = 1 , particles
            write(10,fmt=*)"Xe",r(p,1),r(p,2),r(p,3)
        end do
    close(10)
end subroutine

subroutine kinetic(v,particles,kin,mom)
    implicit none
  integer :: particles , p
  double precision :: v(particles, 3) , kin , mom
  kin = 0.d0
  mom = 0.d0
  do p = 1, particles
      kin = kin + 0.5d0 * ((v(p,1)**2.d0+v(p,2)**2.d0+v(p,3)**2.d0))
      mom = mom + (v(p,1)+v(p,2)+v(p,3))
  end do
end subroutine

subroutine temperature(v,particles,T_instant)
    implicit none
    integer:: particles
    double precision :: T_instant, kin , mom, v(particles,3)
    call kinetic(v,particles,kin,mom)
    T_instant = 2.d0 * kin/(dble(particles)*3.d0)
end subroutine
