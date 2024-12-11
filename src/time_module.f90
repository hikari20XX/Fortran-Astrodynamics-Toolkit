module time_module
    implicit none
    private

    ! Define working precision here if not provided by kind_module:
    integer, parameter :: wp = selected_real_kind(p=15, r=307)

    ! Define conversions locally:
    real(wp), parameter :: day2sec = 86400.0_wp
    real(wp), parameter :: sec2day = 1.0_wp/day2sec

    ! julian date of J2000 epoch
    real(wp), parameter :: jd_j2000 = 2451545.0_wp

    ! Make jd_j2000 public
    public :: jd_j2000

    ! Public interfaces and routines
    public :: julian_day
    public :: julian_date
    public :: julian_date_intsec
    public :: et_to_jd
    public :: jd_to_et
    public :: jd_to_mjd
    public :: mjd_to_jd
    public :: julian_date_to_calendar_date
    public :: time_module_test

    interface julian_date
        module procedure :: julian_date_realsec, julian_date_intsec
    end interface

    interface julian_date_to_calendar_date
        module procedure :: calendar_date_realsec
    end interface

contains

    pure function et_to_jd(et) result(jd)
        ! Convert ephemeris time (sec from J2000) to Julian date
        real(wp), intent(in) :: et
        real(wp) :: jd
        jd = jd_j2000 + et*sec2day
    end function et_to_jd

    pure function jd_to_et(jd) result(et)
        ! Convert Julian date to ephemeris time (sec from J2000)
        real(wp), intent(in) :: jd
        real(wp) :: et
        et = (jd - jd_j2000)*day2sec
    end function jd_to_et

    pure function jd_to_mjd(jd) result(mjd)
        ! Convert Julian date to Modified Julian date
        real(wp), intent(in) :: jd
        real(wp) :: mjd
        mjd = jd - 2400000.5_wp
    end function jd_to_mjd

    pure function mjd_to_jd(mjd) result(jd)
        ! Convert Modified Julian date to Julian date
        real(wp), intent(in) :: mjd
        real(wp) :: jd
        jd = mjd + 2400000.5_wp
    end function mjd_to_jd

    pure integer function julian_day(y,m,d)
        ! Returns the Julian day number for given Y,M,D
        integer,intent(in):: y,m,d
        julian_day = d-32075+1461*(y+4800+(m-14)/12)/4+367*(m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4
    end function julian_day

    pure function julian_date_intsec(y,m,d,hour,minute,second) result(julian_date)
        ! Returns JD for given Y,M,D,H,M,S (integer seconds)
        integer,intent(in):: y,m,d,hour,minute,second
        real(wp):: julian_date
        julian_date = julian_date_realsec(y,m,d,hour,minute,real(second,wp))
    end function julian_date_intsec

    pure function julian_date_realsec(y,m,d,hour,minute,second) result(julian_date)
        ! Returns JD for given Y,M,D,H,M,S (real seconds)
        integer,intent(in):: y,m,d,hour,minute
        real(wp),intent(in):: second
        real(wp):: julian_date
        integer :: julian_day_number
        julian_day_number = julian_day(y,m,d)
        julian_date = real(julian_day_number,wp)+(hour-12.0_wp)/24.0_wp+minute/1440.0_wp+second/86400.0_wp
    end function julian_date_realsec

    pure subroutine calendar_date_realsec(julian_date,year,month,day,hrs,min,sec)
        ! Convert JD to calendar date/time
        real(wp),intent(in):: julian_date
        integer,intent(out):: year,month,day,hrs,min
        real(wp),intent(out):: sec
        integer :: i,j,k,l,n,jd
        real(wp) :: frac_day

        jd = int(julian_date)
        l = jd+68569
        n = 4*l/146097
        l = l-(146097*n+3)/4
        i = 4000*(l+1)/1461001
        l = l-1461*i/4+31
        j = 80*l/2447
        k = l-2447*j/80
        l = j/11
        j = j+2-12*l
        i = 100*(n-49)+i+l
        year=i
        month=j
        day=k
        frac_day = julian_date-real(jd,wp)+0.5_wp
        hrs=int(frac_day*24.0_wp)
        min=int((frac_day - hrs/24.0_wp)*1440.0_wp)
        sec=(frac_day - hrs/24.0_wp - min/1440.0_wp)*86400.0_wp
        if(sec==60.0_wp)then
            sec=0.0_wp
            min=min+1
        end if
        if(min==60)then
            min=0
            hrs=hrs+1
        end if
    end subroutine calendar_date_realsec

    subroutine time_module_test()
        real(wp) :: jd, sec
        integer :: year,month,day,hrs,min
        write(*,*) ''
        write(*,*) '---------------'
        write(*,*) ' time_module_test'
        write(*,*) '---------------'
        write(*,*) ''
        jd = julian_date(2000,1,1,12,0,0)
        call calendar_date_realsec(jd,year,month,day,hrs,min,sec)
        write(*,*) 'jd    ', jd
        write(*,*) 'year  ', year
        write(*,*) 'month ', month
        write(*,*) 'day   ', day
        write(*,*) 'hrs   ', hrs
        write(*,*) 'min   ', min
        write(*,*) 'sec   ', sec
        if (year/=2000)  error stop 'error: incorrect year'
        if (month/=1)    error stop 'error: incorrect month'
        if (day/=1)      error stop 'error: incorrect day'
        if (hrs/=12)     error stop 'error: incorrect hrs'
        if (min/=0)      error stop 'error: incorrect min'
        if (sec/=0.0_wp) error stop 'error: incorrect sec'
        write(*,*) 'PASSED'
    end subroutine time_module_test

end module time_module
