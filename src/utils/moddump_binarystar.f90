!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  Input is a relaxed star, output is two relaxed stars in binary orbit
!
!  REFERENCES: None
!
!  OWNER: Terrence Tricco
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, dim, initial_params, part, prompting, units
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use units,          only: set_units,udist,unit_velocity
 use prompting,      only: prompt
 use centreofmass,   only: reset_centreofmass
 use initial_params, only: get_conserv
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 integer :: opt
 real :: sep,mtot,velocity,corot_vel
 real :: x1com(3), v1com(3), x2com(3), v2com(3)

 print *, 'Running moddump_binarystar: set up binary star systems in close contact'
 print *, ''
 print *, 'Options:'
 print *, '   1) Duplicate a relaxed star'
 print *, '   2) Adjust separation of existing binary'

 opt = 1
 call prompt('Choice',opt, 1, 2)

 if (opt  /=  1 .and. opt  /=  2) then
    print *, 'Incorrect option selected. Doing nothing.'
    return
 endif

 sep = 10.0
 print *, ''
 print *, 'Distance unit is: ', udist
 call prompt('Enter radial separation between stars (in code unit)', sep, 0.)
 print *, ''

 ! duplicate star if chosen
 if (opt == 1) then
    call duplicate_star(npart, npartoftype, massoftype, xyzh, vxyzu)
 endif


 ! add a uniform low density background fluid
! if (opt == 3) then
!    call add_background(npart, npartoftype, massoftype, xyzh, vxyzu)
! endif


 mtot = npart*massoftype(igas)
 velocity  = 0.5 * sqrt(1.0 * mtot) / sqrt(sep) ! in code units
 corot_vel = 2.0 * velocity / sep 

 ! find the centre of mass position and velocity for each star
 call calc_coms(npart,npartoftype,massoftype,xyzh,vxyzu,x1com,v1com,x2com,v2com)
 ! adjust seperation of binary
 call adjust_sep(npart, npartoftype, massoftype, xyzh, vxyzu, sep, x1com, v1com, x2com, v2com)

 ! set orbital velocity
 call set_velocity(npart, npartoftype, massoftype, xyzh, vxyzu, velocity)
 !call set_corotate_velocity(corot_vel)

 ! reset centre of mass of the binary system
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 get_conserv = 1.

end subroutine modify_dump

subroutine duplicate_star(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas,temperature
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 use dim,          only: store_temperature
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real :: sep,mtot,velocity

 npart = npartoftype(igas)

 sep = 10.0

 ! duplicate relaxed star
 do i = npart+1, 2*npart
    ! place star a distance rad away
    xyzh(1,i) = xyzh(1,i-npart) + sep
    xyzh(2,i) = xyzh(2,i-npart)
    xyzh(3,i) = xyzh(3,i-npart)
    xyzh(4,i) = xyzh(4,i-npart)
    vxyzu(1,i) = vxyzu(1,i-npart)
    vxyzu(2,i) = vxyzu(2,i-npart)
    vxyzu(3,i) = vxyzu(3,i-npart)
    vxyzu(4,i) = vxyzu(4,i-npart)
    if (store_temperature) then
       temperature(i) = temperature(i-npart)
    endif
    call set_particle_type(i,igas)
 enddo

 npart = 2 * npart
 npartoftype(igas) = npart

end subroutine duplicate_star

subroutine calc_coms(npart,npartoftype,massoftype,xyzh,vxyzu,x1com,v1com,x2com,v2com)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas,iamtype,iphase,maxphase,maxp
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(out)   :: x1com(:),v1com(:),x2com(:),v2com(:)
 integer :: i, itype
 real    :: xi, yi, zi, vxi, vyi, vzi
 real    :: totmass, pmassi, dm
 ! calc centre of mass of each star to form the reference points to adjust the position of the second star

 ! first star
 x1com = 0.
 v1com = 0.
 totmass = 0.
 do i = 1, npart/2
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    if (maxphase == maxp) then
       itype = iamtype(iphase(i))
       if (itype > 0) then
          pmassi = massoftype(itype)
       else
          pmassi = massoftype(igas)
       endif
    else
       pmassi = massoftype(igas)
    endif

    totmass = totmass + pmassi
    x1com(1) = x1com(1) + pmassi * xi
    x1com(2) = x1com(2) + pmassi * yi
    x1com(3) = x1com(3) + pmassi * zi
    v1com(1) = v1com(1) + pmassi * vxi
    v1com(2) = v1com(2) + pmassi * vyi
    v1com(3) = v1com(3) + pmassi * vzi
 enddo

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 x1com = dm * x1com
 v1com = dm * v1com

 ! second star
 x2com = 0.
 v2com = 0.
 totmass = 0.
 do i = npart/2 + 1, npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)
    if (maxphase == maxp) then
       itype = iamtype(iphase(i))
       if (itype > 0) then
          pmassi = massoftype(itype)
       else
          pmassi = massoftype(igas)
       endif
    else
       pmassi = massoftype(igas)
    endif

    totmass = totmass + pmassi
    x2com(1) = x2com(1) + pmassi * xi
    x2com(2) = x2com(2) + pmassi * yi
    x2com(3) = x2com(3) + pmassi * zi
    v2com(1) = v2com(1) + pmassi * vxi
    v2com(2) = v2com(2) + pmassi * vyi
    v2com(3) = v2com(3) + pmassi * vzi
 enddo

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 x2com = dm * x2com
 v2com = dm * v2com

end subroutine calc_coms

subroutine adjust_sep(npart,npartoftype,massoftype,xyzh,vxyzu,sep,x1com,v1com,x2com,v2com)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: x1com(:),v1com(:),x2com(:),v2com(:)
 real,    intent(in)    :: sep
 integer :: i

 ! now we now the centre point of each star, we can set star 1 to origin, star 2 sep away on x axis, then reset com
 do i = 1, npart/2
    xyzh(1,i) = xyzh(1,i) - x1com(1)
    xyzh(2,i) = xyzh(2,i) - x1com(2)
    xyzh(3,i) = xyzh(3,i) - x1com(3)
    vxyzu(1,i) = vxyzu(1,i) - v1com(1)
    vxyzu(2,i) = vxyzu(2,i) - v1com(2)
    vxyzu(3,i) = vxyzu(3,i) - v1com(3)
 enddo

 do i = npart/2 + 1, npart
    xyzh(1,i) = xyzh(1,i) - x2com(1) + sep
    xyzh(2,i) = xyzh(2,i) - x2com(2)
    xyzh(3,i) = xyzh(3,i) - x2com(3)
    vxyzu(1,i) = vxyzu(1,i) - v2com(1)
    vxyzu(2,i) = vxyzu(2,i) - v2com(2)
    vxyzu(3,i) = vxyzu(3,i) - v2com(3)
 enddo

end subroutine adjust_sep


subroutine set_corotate_velocity(corot_vel)
 use options,        only:iexternalforce
 use externalforces, only: omega_corotate,iext_corotate
 real,    intent(in)    :: corot_vel
 integer :: i

 !turns on corotation
 iexternalforce = iext_corotate
 omega_corotate = corot_vel

 print "(/,a,es18.10,/)", ' The angular velocity in the corotating frame is: ', omega_corotate
end subroutine


subroutine set_velocity(npart,npartoftype,massoftype,xyzh,vxyzu,velocity)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass,igas,set_particle_type,igas
 use units,        only: set_units,udist,unit_velocity
 use prompting,    only: prompt
 use centreofmass, only: reset_centreofmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: velocity
 integer :: i
 real :: mtot

 print *, "Adding a bulk velocity |V| = ", velocity, "( = ", (velocity*unit_velocity), &
                  " physical units) to set stars in mutual orbit"
 print *, ''
 ! Adjust bulk velocity of relaxed star towards second star
 do i = 1, npart/2
    vxyzu(2,i) = vxyzu(2,i) + velocity
 enddo

 do i = npart/2 + 1, npart
    vxyzu(2,i) = vxyzu(2,i) - velocity
 enddo

end subroutine set_velocity


end module moddump

