program q2_numeric_implicit
    implicit none
    
    ! --- Variable Declarations ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: m  = 6  ! Number of isotopes (U-238 -> ... -> Pb-206)
    
    ! Isotope Data
    character(len=15), parameter :: names(m) = [ &
        'U-238          ', &
        'U-234          ', &
        'Th-230         ', &
        'Ra-226         ', &
        'Pb-210         ', &
        'Pb-206 (stable)']
        
    real(dp), parameter :: A_molar(m) = [238.0_dp, 234.0_dp, 230.0_dp, 226.0_dp, 210.0_dp, 206.0_dp]
    real(dp), parameter :: t12(m)     = [4.5e9_dp, 2.25e5_dp, 8.0e4_dp, 1.6e3_dp, 22.0_dp, huge(1.0_dp)]
    
    real(dp) :: lam(m), ln2
    integer  :: i, j
    
    ! Constants
    real(dp), parameter :: NA       = 6.02214076e23_dp
    real(dp), parameter :: G_PER_KG = 1000.0_dp
    real(dp), parameter :: T_END    = 4.503e9_dp
    
    ! Initial Condition (Mass of U-238 today)
    real(dp), parameter :: M_today_U238_kg = 2.892e15_dp
    
    ! Arrays for calculation
    real(dp) :: N0(m), N_T_num(m), N_T_ana(m)
    real(dp) :: M_T_num(m), M_T_ana(m)
    real(dp) :: N_today_U238, N0_parent
    
    ! Time steps to test
    real(dp), parameter :: dts(3) = [1.0e6_dp, 1.0e4_dp, 1.0e2_dp]
    integer :: nsteps
    
    ! --- Initialization ---
    ln2 = log(2.0_dp)
    do i = 1, m
        if (t12(i) < huge(1.0_dp)) then
            lam(i) = ln2 / t12(i)
        else
            lam(i) = 0.0_dp
        end if
    end do
    
    ! 1. Calculate Initial Parent (N0) at t=0 based on Today's U-238
    N_today_U238 = masskg_to_atoms(M_today_U238_kg, A_molar(1))
    N0_parent    = N_today_U238 * exp(lam(1) * T_END)
    
    ! Set initial vector (only U-238 exists at t=0)
    N0 = 0.0_dp
    N0(1) = N0_parent
    
    ! Calculate Analytic Solution for Reference
    call bateman_solution(T_END, lam, N0_parent, N_T_ana)
    
    ! Convert analytic results to mass
    do i = 1, m
        M_T_ana(i) = atoms_to_masskg(N_T_ana(i), A_molar(i))
    end do

    print *, '--- IMPLICIT EULER SIMULATION ---'
    
    ! --- Main Loop for Different Time Steps ---
    do j = 1, 3
        print *, ' '
        print '(A, ES10.1, A)', '>>> Simulation with dt = ', dts(j), ' years'
        
        call implicit_euler_forward(N0, dts(j), nsteps, N_T_num)
        
        ! Convert results to mass
        do i = 1, m
            M_T_num(i) = atoms_to_masskg(N_T_num(i), A_molar(i))
        end do
        
        print '(A, I0)', '    Steps: ', nsteps
        print '(A)', '    Isotope            Numeric (kg)       Analytic (kg)      Rel. Err (%)'
        print '(A)', '    ----------------------------------------------------------------'
        
        do i = 1, m
            print '(4X, A15, 2(ES16.6), F12.4)', &
                  names(i), M_T_num(i), M_T_ana(i), relerr_pct(M_T_num(i), M_T_ana(i))
        end do
    end do

contains

    ! --- Helper Functions ---

    pure real(dp) function masskg_to_atoms(mass_kg, M_gpmol) result(N)
        real(dp), intent(in) :: mass_kg, M_gpmol
        N = (mass_kg * G_PER_KG / M_gpmol) * NA
    end function

    pure real(dp) function atoms_to_masskg(N, M_gpmol) result(mkg)
        real(dp), intent(in) :: N, M_gpmol
        mkg = (N / NA) * M_gpmol / G_PER_KG
    end function

    pure real(dp) function relerr_pct(num, ref) result(p)
        real(dp), intent(in) :: num, ref
        if (ref == 0.0_dp) then
            p = 0.0_dp
        else
            p = 100.0_dp * abs(num - ref) / ref
        end if
    end function

    ! --- Solvers ---

    subroutine bateman_solution(Time, lambdas, N_initial, N_out)
        real(dp), intent(in)  :: Time, N_initial
        real(dp), intent(in)  :: lambdas(m)
        real(dp), intent(out) :: N_out(m)
        real(dp) :: num, total, denom
        integer  :: r, j2, k2

        N_out = 0.0_dp
        N_out(1) = N_initial * exp(-lambdas(1) * Time)
        
        do k2 = 2, m-1
            num = 1.0_dp
            do r = 1, k2-1
                num = num * lambdas(r)
            end do
            
            total = 0.0_dp
            do j2 = 1, k2
                denom = 1.0_dp
                do r = 1, k2
                    if (r /= j2) denom = denom * (lambdas(r) - lambdas(j2))
                end do
                total = total + exp(-lambdas(j2) * Time) / denom
            end do
            N_out(k2) = N_initial * num * total
        end do
        
        ! Conservation for stable element
        N_out(m) = N_initial
        do k2 = 1, m-1
            N_out(m) = N_out(m) - N_out(k2)
        end do
    end subroutine

    subroutine implicit_euler_forward(N_in, dt, steps_out, N_out)
        real(dp), intent(in)  :: N_in(m), dt
        integer,  intent(out) :: steps_out
        real(dp), intent(out) :: N_out(m)
        integer  :: step, k
        real(dp) :: denominator(m)
        
        steps_out = ceiling(T_END / dt)
        N_out = N_in
        
        ! Pre-calculate denominators (1 + lambda*dt) for efficiency
        denominator = 1.0_dp
        do k = 1, m-1
            denominator(k) = 1.0_dp / (1.0_dp + lam(k)*dt)
        end do
        
        do step = 1, steps_out
            ! Parent: N1 = N1_prev / (1 + lam1*dt)
            N_out(1) = N_out(1) * denominator(1)
            
            ! Daughters: Nk = (Nk_prev + Production) / (1 + lamk*dt)
            do k = 2, m-1
                N_out(k) = (N_out(k) + lam(k-1)*dt*N_out(k-1)) * denominator(k)
            end do
            
            ! Stable Product: N_stable += Production
            N_out(m) = N_out(m) + lam(m-1)*dt*N_out(m-1)
        end do
    end subroutine

end program q2_numeric_implicit
