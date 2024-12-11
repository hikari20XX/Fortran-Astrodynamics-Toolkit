module time_module
    use kind_module
    implicit none
    private

    real(wp), parameter, public :: jd_j2000 = 2451545.0_wp

    interface julian_date
        module procedure julian_date_realsec
    end interface

    public :: julian_day
    public :: julian_date_realsec
    public :: et_to_jd
    public :: jd_to_et
    public :: jd_to_mjd
    public :: mjd_to_jd
    public :: julian_date_to_calendar_date

contains

    pure function et_to_jd(et) result(jd)
        use conversion_module, only: sec2day
        implicit none
        real(wp), intent(in) :: et
        real(wp) :: jd

        jd = jd_j2000 + et * sec2day
    end function et_to_jd

    pure function jd_to_et(jd) result(et)
        use conversion_module, only: day2sec
        implicit none
        real(wp), intent(in) :: jd
        real(wp) :: et

        et = (jd - jd_j2000)*day2sec
    end function jd_to_et

    pure function jd_to_mjd(jd) result(mjd)
        implicit none
        real(wp), intent(in) :: jd
        real(wp) :: mjd
        mjd = jd - 2400000.5_wp
    end function jd_to_mjd

    pure function mjd_to_jd(mjd) result(jd)
        implicit none
        real(wp), intent(in) :: mjd
        real(wp) :: jd
        jd = mjd + 2400000.5_wp
    end function mjd_to_jd

    pure integer function julian_day(y,m,d)
        implicit none
        integer,intent(in):: y,m,d
        julian_day = d-32075+1461*(y+4800+(m-14)/12)/4+367*(m-2-(m-14)/12*12)/12 - 3*((y+4900+(m-14)/12)/100)/4
    end function julian_day

    pure function julian_date_realsec(y,m,d,hour,minute,second) result(julian_date)
        implicit none
        integer,intent(in):: y,m,d,hour,minute
        real(wp),intent(in):: second
        real(wp):: julian_date
        integer :: jdn

        jdn = julian_day(y,m,d)
        julian_date = real(jdn,wp) + (hour-12.0_wp)/24.0_wp + minute/1440.0_wp + second/86400.0_wp
    end function julian_date_realsec

    pure subroutine julian_date_to_calendar_date(jd,year,month,day,hrs,min,sec)
        implicit none
        real(wp),intent(in):: jd
        integer,intent(out):: year,month,day,hrs,min
        real(wp),intent(out):: sec

        integer :: i,j,k,l,n,jdi
        real(wp):: frac_day
        jdi = int(jd)
        l = jdi+68569
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

        frac_day = jd - real(jdi,wp) + 0.5_wp
        hrs=int(frac_day*24.0_wp)
        min=int((frac_day - hrs/24.0_wp)*1440.0_wp)
        sec=(frac_day - hrs/24.0_wp - min/1440.0_wp)*86400.0_wp

        if (sec==60.0_wp) then
           sec=0.0_wp
           min=min+1
        end if

        if (min==60) then
           min=0
           hrs=hrs+1
        end if
    end subroutine julian_date_to_calendar_date

end module time_module
