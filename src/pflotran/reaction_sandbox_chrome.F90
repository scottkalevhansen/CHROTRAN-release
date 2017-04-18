module Reaction_Sandbox_Chrome_class

! Sandbox reaction for Cr(VI) reduction using bio-reduction

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_Cr_bio_type
    character(len=MAXWORDLENGTH) :: food_name   ! This the food in the reaction
    character(len=MAXWORDLENGTH) :: Cr_name     ! This is for Cr(VI) in the reaction
    character(len=MAXWORDLENGTH) :: biomass_name ! This is for the biomass in the reaction
    PetscInt :: food_id 
    PetscInt :: Cr_id 
    PetscInt :: biomass_id
    PetscReal :: food_decay_constant_1
    PetscReal :: food_decay_constant_2
    PetscReal :: Cr_decay_constant
    PetscReal :: Cr_inhibition_constant
    PetscReal :: biomass_growth_constant
    PetscReal :: biomass_decay_constant 
    PetscReal :: monod_food_constant
    PetscReal :: inhibition_biomass_constant
    PetscReal :: biomass_min    ! Minimum background concentration of the biomass
  contains
    procedure, public :: ReadInput => ChromeRead
    procedure, public :: Setup => ChromeSetup
    procedure, public :: Evaluate => ChromeReact
    procedure, public :: Destroy => ChromeDestroy
  end type reaction_sandbox_Cr_bio_type

  public :: ChromeCreate

contains

! ************************************************************************** !

function ChromeCreate()
  ! 
  ! Allocates Cr(VI) bio-reduction object
  ! 
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  ! 

  implicit none
  
  class(reaction_sandbox_Cr_bio_type), pointer :: ChromeCreate

  allocate(ChromeCreate)
  ChromeCreate%food_name = ''
  ChromeCreate%Cr_name = ''
  ChromeCreate%biomass_name = ''
  ChromeCreate%food_id = 0
  ChromeCreate%Cr_id = 0
  ChromeCreate%biomass_id = 0
  ChromeCreate%food_decay_constant_1 = 0.d0
  ChromeCreate%food_decay_constant_2 = 0.d0
  ChromeCreate%Cr_decay_constant = 0.d0
  ChromeCreate%Cr_inhibition_constant = 0.d0
  ChromeCreate%biomass_decay_constant = 0.d0
  ChromeCreate%biomass_growth_constant = 0.d0
  ChromeCreate%monod_food_constant = 0.d0
  ChromeCreate%inhibition_biomass_constant = 0.d0
  ChromeCreate%biomass_min = 0.d0
  nullify(ChromeCreate%next)  
      
end function ChromeCreate

! ************************************************************************** !

subroutine ChromeRead(this,input,option)
  ! 
  ! Reads input deck for Cr(VI) bio-reduction parameters 
  ! 
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  
  implicit none
  
  class(reaction_sandbox_Cr_bio_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
    call StringToUpper(word)   

    select case(trim(word))

      ! Chrome Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     CR-BIO
      !       FOOD_NAME
      !       CR_NAME 
      !       BIOMASS_NAME 
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

      case('FOOD_NAME')
        call InputReadWord(input,option,this%food_name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'food_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')    
      case('CR_NAME')
        call InputReadWord(input,option,this%Cr_name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'Cr_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')    
      case('BIOMASS_NAME')
        call InputReadWord(input,option,this%biomass_name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'biomass_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')    
      case('FOOD_DECAY_CONSTANT_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%food_decay_constant_1)
        call InputErrorMsg(input,option,'food_decay_constant_1', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('FOOD_DECAY_CONSTANT_2')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%food_decay_constant_2)
        call InputErrorMsg(input,option,'food_decay_constant_2', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('CR_DECAY_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%Cr_decay_constant)
        call InputErrorMsg(input,option,'Cr_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('CR_INHIBITION_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%Cr_inhibition_constant)
        call InputErrorMsg(input,option,'Cr_inhibition_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('BIOMASS_GROWTH_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%biomass_growth_constant)
        call InputErrorMsg(input,option,'biomass_growth_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
     case('BIOMASS_DECAY_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%biomass_decay_constant)
        call InputErrorMsg(input,option,'biomass_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
     case('MONOD_FOOD_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%monod_food_constant)
        call InputErrorMsg(input,option,'monod_food_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('INHIBITION_BIOMASS_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_biomass_constant)
        call InputErrorMsg(input,option,'inhibition_biomass_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case('BIOMASS_BACKGROUND_CONCENTRATION')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%biomass_min )
        call InputErrorMsg(input,option,'biomass_background_concentration', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO')
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,CR-BIO',option)
    end select
  enddo
  
end subroutine ChromeRead

! ************************************************************************** !

subroutine ChromeSetup(this,reaction,option)
  ! 
  ! Sets up the Cr(VI) bio-reduction reaction 
  ! 
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_Cr_bio_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  this%food_id = &
    GetPrimarySpeciesIDFromName(this%food_name, &
                                reaction,option)

  this%Cr_id = &
    GetPrimarySpeciesIDFromName(this%Cr_name, &
                                reaction,option)
 
  this%biomass_id = &
    GetImmobileSpeciesIDFromName(this%biomass_name, &
                                 reaction%immobile,option)

end subroutine ChromeSetup

! ************************************************************************** !

subroutine ChromeReact(this,Residual,Jacobian,compute_derivative, &
                       rt_auxvar,global_auxvar,material_auxvar, &
                       reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_Cr_bio_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: mu, dmu_dF, dmu_dB
  PetscInt :: idof_food, idof_biomass, idof_Cr
  
  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives 
  !   or a numerical approximation unless always running with 
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of 
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error 
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in Chrome'
  !   call printErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.  
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations 
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3] 
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3] 
  !     global_auxvar%sat(iphase) - saturation 
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]                     
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.
  
  ! Unit of the residual must be in moles/second
  ! global_auxvar%sat(iphase) = saturation of cell
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! alway subtract contribution from residual

  idof_food = this%food_id
  idof_Cr = this%Cr_id
  idof_biomass = reaction%offset_immobile + this%biomass_id

#if 0
  print *, 'food_decay_constant_1:', this%food_decay_constant_1 
  print *, 'food_decay_constant_2:', this%food_decay_constant_2 
  print *, 'Cr_decay_constant:', this%Cr_decay_constant
  print *, 'biomass_growth_constant:', this%biomass_growth_constant
  print *, 'biomass_decay_constant:', this%biomass_decay_constant
  print *, 'monod_food_constant:', this%monod_food_constant
  print *, 'inhibition_biomass_constant:', this%inhibition_biomass_constant

  print *, 'immobile:', rt_auxvar%immobile(this%biomass_id)
  print *, 'food:', rt_auxvar%total(idof_food,iphase)
  print *, 'Cr:', rt_auxvar%total(idof_Cr,iphase)
#endif 
 
  mu = this%biomass_growth_constant*rt_auxvar%immobile(this%biomass_id)* &
       (rt_auxvar%total(idof_food,iphase)/ &
       (rt_auxvar%total(idof_food,iphase) + this%monod_food_constant))* &
       (this%inhibition_biomass_constant/ &
       (rt_auxvar%immobile(this%biomass_id) + this%inhibition_biomass_constant))

#if 0  
  print *, 'mu:', mu
#endif 

  Residual(idof_food) = Residual(idof_food) + &
                           this%food_decay_constant_1* &
                           mu*material_auxvar%volume + &
                           this%food_decay_constant_2* &                          ! 1/s 
                           rt_auxvar%immobile(this%biomass_id)* &                 ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk 
  Residual(idof_Cr) = Residual(idof_Cr) + &
                         this%Cr_decay_constant/ &                                ! L/mol/s
                         (this%Cr_inhibition_constant + &
                          rt_auxvar%total(idof_Cr,iphase))* &
                         rt_auxvar%total(idof_Cr,iphase)* &                       ! mol/L
                         rt_auxvar%immobile(this%biomass_id)* &                   ! mol/m3 bulk
                         material_auxvar%volume                                   ! m3 bulk
 
  Residual(idof_biomass) = Residual(idof_biomass) - &
                           mu*material_auxvar%volume + &                          ! mol/m3 bulk/s * m3 bulk
                           this%biomass_decay_constant* &                         ! 1/s
                           (rt_auxvar%immobile(this%biomass_id) - &
                            this%biomass_min)* &                                  ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk

#if 0
  print *, 'Residual(idof_food):', Residual(idof_food)
  print *, 'Residual(idof_Cr):', Residual(idof_Cr)
  print *, 'Residual(idof_biomass):', Residual(idof_biomass)
#endif 
  
  if (compute_derivative) then


    dmu_dF = this%biomass_growth_constant* &
             rt_auxvar%immobile(this%biomass_id)* &
             (this%monod_food_constant/ &
              (rt_auxvar%total(idof_food,iphase) + &
               this%monod_food_constant)**2)* &
             (this%inhibition_biomass_constant/ &
              (rt_auxvar%immobile(this%biomass_id) + &
               this%inhibition_biomass_constant))
             

    dmu_dB = this%biomass_growth_constant* &
             (rt_auxvar%total(idof_food,iphase)/ &
              (rt_auxvar%total(idof_food,iphase) + &
               this%monod_food_constant))* &
             (this%inhibition_biomass_constant**2/ &
              (rt_auxvar%immobile(this%biomass_id) + &
               this%inhibition_biomass_constant)**2)


    Jacobian(idof_food,idof_food) =  &
    Jacobian(idof_food,idof_food) + &
      this%food_decay_constant_1* &                            ! unitless 
      dmu_dF* &                                                ! L water / m3 bulk / s
      material_auxvar%volume* &                                ! m3 bulk
      rt_auxvar%aqueous%dtotal(idof_food,idof_food,iphase)     ! kg water / L water

    
    ! Note that this part of the Jacobian needs to be in m3 bulk / s 
    Jacobian(idof_food,idof_biomass) = &
    Jacobian(idof_food,idof_biomass) - &
      (-this%food_decay_constant_2 - &
       this%food_decay_constant_1*dmu_dB)* &   ! 1/s
      material_auxvar%volume                                                   ! m3 bulk

    Jacobian(idof_Cr,idof_Cr) = &
    Jacobian(idof_Cr,idof_Cr) + &
      this%Cr_decay_constant/ &                                                ! L/mol/s
     (this%Cr_inhibition_constant + &
      rt_auxvar%total(idof_Cr,iphase))* &
     (1.d0 - rt_auxvar%total(idof_Cr,iphase)/ &
      (rt_auxvar%total(idof_Cr,iphase) + this%Cr_inhibition_constant))* &
      rt_auxvar%immobile(this%biomass_id)* &                                   ! mol/m3 bulk
      material_auxvar%volume* &                                                ! m3 bulk
      rt_auxvar%aqueous%dtotal(idof_Cr,idof_Cr,iphase)                         ! kg water/L water

    ! Note that this part of the Jacobian needs to be in m3 bulk / s 
    Jacobian(idof_Cr,idof_biomass) = &
    Jacobian(idof_Cr,idof_biomass) + &
      this%Cr_decay_constant/ &                                                ! L/mol/s
     (this%Cr_inhibition_constant + &
      rt_auxvar%total(idof_Cr,iphase))* &
      rt_auxvar%total(idof_Cr,iphase)* &                                       ! mol/L
      material_auxvar%volume                                                   ! m3 bulk

    Jacobian(idof_biomass,idof_food) = &
    Jacobian(idof_biomass,idof_food) - &
      dmu_dF* &                                                                ! L water/m3 bulk/s
      material_auxvar%volume* &                                                ! m3 bulk
      rt_auxvar%aqueous%dtotal(idof_food,idof_food,iphase)                     ! kg water/L water


    Jacobian(idof_biomass,idof_biomass) = & 
    Jacobian(idof_biomass,idof_biomass) - & 
      (dmu_dB - this%biomass_decay_constant)* &                             ! 1/s 
      material_auxvar%volume                                                ! m3 bulk

  endif

  
end subroutine ChromeReact

! ************************************************************************** !

subroutine ChromeDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  ! 

  implicit none
  
  class(reaction_sandbox_Cr_bio_type) :: this  

! 12. Add code to deallocate contents of the Cr_bio object

end subroutine ChromeDestroy

end module Reaction_Sandbox_Chrome_class
