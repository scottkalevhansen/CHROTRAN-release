module Reaction_Sandbox_Chrome_Alcohol_Perm_class

! Sandbox reaction for Cr(VI) reduction using bio-reduction with reduced permeability

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_Cr_bio_perm_type
    character(len=MAXWORDLENGTH) :: food_mobile_name   ! This the mobile food in the reaction
    character(len=MAXWORDLENGTH) :: food_immobile_name ! This the immobile food in the reaction
    character(len=MAXWORDLENGTH) :: Cr_name     ! This is for Cr(VI) in the reaction
    character(len=MAXWORDLENGTH) :: biomass_name ! This is for the biomass in the reaction
    character(len=MAXWORDLENGTH) :: alcohol_name ! This is for the alcohol in the reaction
    character(len=MAXWORDLENGTH) :: biocide_name ! This is for the biocide in the reaction
    character(len=MAXWORDLENGTH) :: biomineral_name ! This is for the dummny bio mineral

    PetscInt :: food_immobile_id  ! Immobile food
    PetscInt :: food_mobile_id    ! Mobile food
    PetscInt :: Cr_id
    PetscInt :: biomass_id
    PetscInt :: biovolume_id
    PetscInt :: alcohol_id
    PetscInt :: biocide_id
    PetscInt :: biomineral_id

    ! Decay and inhibition parameters in our sophisticated model
    PetscReal :: food_decay_constant_1
    PetscReal :: food_decay_constant_2
    PetscReal :: Cr_decay_constant
    PetscReal :: Cr_inhibition_constant
    PetscReal :: biomass_growth_constant
    PetscReal :: biomass_decay_constant
    PetscReal :: monod_food_constant
    PetscReal :: inhibition_biomass_constant
    PetscReal :: biomass_min    ! Minimum background concentration of the biomass
    PetscReal :: Cr_food_direct_constant
    PetscReal :: Cr_food_Cr_stoic_constant
    PetscReal :: Cr_food_food_stoic_constant
    PetscReal :: food_immobilization_constant
    PetscReal :: food_remobilization_constant
    PetscReal :: inhibition_biomass_exponent
    PetscReal :: inhibition_alcohol_constant
    PetscReal :: biocide_decay_constant
    PetscReal :: biomass_biocide_decay_constant
    PetscReal :: biomass_density
    
  contains
    procedure, public :: ReadInput => ChromeAlcoholPermRead
    procedure, public :: Setup => ChromeAlcoholPermSetup
    procedure, public :: Evaluate => ChromeAlcoholPermReact
    procedure, public :: UpdateKineticState => ChromeAlcoholPermKineticState
    procedure, public :: Destroy => ChromeAlcoholPermDestroy
  end type reaction_sandbox_Cr_bio_perm_type

  public :: ChromeAlcoholPermCreate

contains

! ************************************************************************** !

function ChromeAlcoholPermCreate()
  !
  ! Allocates Cr(VI) bio-reduction object
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  implicit none

  class(reaction_sandbox_Cr_bio_perm_type), pointer :: ChromeAlcoholPermCreate

  allocate(ChromeAlcoholPermCreate)
  ChromeAlcoholPermCreate%food_mobile_name = ''
  ChromeAlcoholPermCreate%food_immobile_name = ''
  ChromeAlcoholPermCreate%Cr_name = ''
  ChromeAlcoholPermCreate%biomass_name = ''
  ChromeAlcoholPermCreate%alcohol_name = ''
  ChromeAlcoholPermCreate%biocide_name = ''
  ChromeAlcoholPermCreate%biomineral_name = ''
! 
  ChromeAlcoholPermCreate%food_mobile_id = 0
  ChromeAlcoholPermCreate%food_immobile_id = 0
  ChromeAlcoholPermCreate%Cr_id = 0
  ChromeAlcoholPermCreate%biomass_id = 0
  ChromeAlcoholPermCreate%alcohol_id = 0
  ChromeAlcoholPermCreate%biocide_id = 0
  ChromeAlcoholPermCreate%biomineral_id = 0
  
  ChromeAlcoholPermCreate%food_decay_constant_1 = 0.d0
  ChromeAlcoholPermCreate%food_decay_constant_2 = 0.d0
  ChromeAlcoholPermCreate%Cr_decay_constant = 0.d0
  ChromeAlcoholPermCreate%Cr_inhibition_constant = 0.d0
  ChromeAlcoholPermCreate%biomass_decay_constant = 0.d0
  ChromeAlcoholPermCreate%biomass_growth_constant = 0.d0
  ChromeAlcoholPermCreate%monod_food_constant = 0.d0
  ChromeAlcoholPermCreate%inhibition_biomass_constant = 0.d0
  ChromeAlcoholPermCreate%biomass_min = 0.d0
  ChromeAlcoholPermCreate%Cr_food_direct_constant = 0.d0
  ChromeAlcoholPermCreate%Cr_food_Cr_stoic_constant = 1.d0
  ChromeAlcoholPermCreate%Cr_food_food_stoic_constant = 1.d0
  ChromeAlcoholPermCreate%food_immobilization_constant = 0.d0
  ChromeAlcoholPermCreate%food_remobilization_constant = 0.d0
  ChromeAlcoholPermCreate%inhibition_biomass_exponent = 0.d0
  ChromeAlcoholPermCreate%inhibition_alcohol_constant = 0.d0
  ChromeAlcoholPermCreate%biocide_decay_constant = 0.d0
  ChromeAlcoholPermCreate%biomass_biocide_decay_constant = 0.d0
  ChromeAlcoholPermCreate%biomass_density = 0.d0
  nullify(ChromeAlcoholPermCreate%next)

end function ChromeAlcoholPermCreate

! ************************************************************************** !

subroutine ChromeAlcoholPermRead(this,input,option)
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

  class(reaction_sandbox_Cr_bio_perm_type) :: this
  type(input_type), pointer  :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL-PERM')
    call StringToUpper(word)

    select case(trim(word))

      ! ChromeAlcohol Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     CR-BIO-ALCOHOL
      !       FOOD_NAME
      !       CR_NAME
      !       BIOMASS_NAME
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

      case('FOOD_MOBILE_NAME')
        call InputReadWord(input,option,this%food_mobile_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'food_mobile_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('FOOD_IMMOBILE_NAME')
        call InputReadWord(input,option,this%food_immobile_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'food_immobile_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_NAME')
        call InputReadWord(input,option,this%Cr_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'Cr_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOMASS_NAME')
        call InputReadWord(input,option,this%biomass_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'biomass_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('ALCOHOL_NAME')
        call InputReadWord(input,option,this%alcohol_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'alcohol_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOCIDE_NAME')
        call InputReadWord(input,option,this%biocide_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'biocide_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')  
      case('BIOMINERAL_NAME')
        call InputReadWord(input,option,this%biomineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'biomineral_name', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')    
      case('FOOD_DECAY_CONSTANT_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%food_decay_constant_1)
        call InputErrorMsg(input,option,'food_decay_constant_1', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('FOOD_DECAY_CONSTANT_2')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%food_decay_constant_2)
        call InputErrorMsg(input,option,'food_decay_constant_2', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_DECAY_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%Cr_decay_constant)
        call InputErrorMsg(input,option,'Cr_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_INHIBITION_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%Cr_inhibition_constant)
        call InputErrorMsg(input,option,'Cr_inhibition_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOMASS_GROWTH_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%biomass_growth_constant)
        call InputErrorMsg(input,option,'biomass_growth_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOMASS_DECAY_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%biomass_decay_constant)
        call InputErrorMsg(input,option,'biomass_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('MONOD_FOOD_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%monod_food_constant)
        call InputErrorMsg(input,option,'monod_food_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('INHIBITION_BIOMASS_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_biomass_constant)
        call InputErrorMsg(input,option,'inhibition_biomass_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOMASS_BACKGROUND_CONCENTRATION')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%biomass_min )
        call InputErrorMsg(input,option,'biomass_background_concentration', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_FOOD_DIRECT_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%Cr_food_direct_constant )
        call InputErrorMsg(input,option,'Cr_food_direct_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_FOOD_CR_STOIC_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%Cr_food_Cr_stoic_constant )
        call InputErrorMsg(input,option,'Cr_food_Cr_stoic_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('CR_FOOD_FOOD_STOIC_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%Cr_food_food_stoic_constant)
        call InputErrorMsg(input,option,'Cr_food_food_stoic_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('FOOD_IMMOBILIZATION_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%food_immobilization_constant)
        call InputErrorMsg(input,option,'food_immobilization_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('FOOD_REMOBILIZATION_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%food_remobilization_constant)
        call InputErrorMsg(input,option,'food_remobilization_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')

      case('INHIBITION_BIOMASS_EXPONENT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%inhibition_biomass_exponent)
        call InputErrorMsg(input,option,'inhibition_biomass_exponent', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('INHIBITION_ALCOHOL_CONSTANT')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%inhibition_alcohol_constant)
        call InputErrorMsg(input,option,'inhibition_alcohol_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
                           
      case('BIOCIDE_DECAY_CONSTANT')
        call InputReadDouble(input,option,this%biocide_decay_constant)
        call InputErrorMsg(input,option,'biocide_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case('BIOMASS_BIOCIDE_DECAY_CONSTANT')
        call InputReadDouble(input,option,this%biomass_biocide_decay_constant)
        call InputErrorMsg(input,option,'biomass_biocide_decay_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')                         
      case('BIOMASS_DENSITY')
        call InputReadDouble(input,option,this%biomass_density)
        call InputErrorMsg(input,option,'biomass_density', &
                           'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL')
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,CR-BIO-ALCOHOL',option)
    end select
  enddo

end subroutine ChromeAlcoholPermRead

! ************************************************************************** !

subroutine ChromeAlcoholPermSetup(this,reaction,option)
  !
  ! Sets up the Cr(VI) bio-reduction reaction
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_Cr_bio_perm_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  this%food_mobile_id = &
    GetPrimarySpeciesIDFromName(this%food_mobile_name, &
                                reaction,option)

  this%Cr_id = &
    GetPrimarySpeciesIDFromName(this%Cr_name, &
                                reaction,option)

  this%alcohol_id = &
    GetPrimarySpeciesIDFromName(this%alcohol_name, &
                                reaction,option)

  this%biocide_id = &
    GetPrimarySpeciesIDFromName(this%biocide_name, &
                                reaction,option)

  this%biomass_id = &
    GetImmobileSpeciesIDFromName(this%biomass_name, &
                                 reaction%immobile,option)
                                 
  this%food_immobile_id = &
    GetImmobileSpeciesIDFromName(this%food_immobile_name, &
                                 reaction%immobile,option)
  this%biomineral_id = &
    GetMineralIDFromName(this%biomineral_name, &
                         reaction%mineral,option)
  
end subroutine ChromeAlcoholPermSetup

! ************************************************************************** !

subroutine ChromeAlcoholPermReact(this,Residual,Jacobian,compute_derivative, &
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

  class(reaction_sandbox_Cr_bio_perm_type) :: this
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
  PetscReal :: mu_B, mu_D
  PetscReal :: dmu_dF, dmu_dB, sum_food
  PetscInt :: idof_food_mobile, idof_food_immobile, idof_biomass, idof_Cr
  PetscInt :: idof_alcohol, idof_biomineral, idof_biocide
  PetscReal :: immobile_to_water_vol
  PetscReal :: immobile_mole_fraction, mobile_mole_fraction
  PetscReal :: biomass_residual_delta

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
  !                      'due to assumptions in ChromeAlcohol'
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
  ! always subtract contribution from residual

  idof_food_mobile = this%food_mobile_id
  idof_Cr = this%Cr_id
  idof_alcohol = this%alcohol_id
  idof_biocide = this%biocide_id
  idof_biomass = reaction%offset_immobile + this%biomass_id
  idof_food_immobile = reaction%offset_immobile + this%food_immobile_id

  immobile_to_water_vol = &
     material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0                   ! L water/ m3 bulk

  sum_food = rt_auxvar%total(idof_food_mobile,iphase) + &
             rt_auxvar%immobile(this%food_immobile_id)/ &
             immobile_to_water_vol                                                ! in mol/L water; Note that food_immobile is divided by porosity*saturation

  mu_B = this%biomass_growth_constant*rt_auxvar%immobile(this%biomass_id)* &      ! mol/m3 bulk/s
        ! F monod term, unitless
        (sum_food/(sum_food + this%monod_food_constant))* &
        ! B monod inhibition term, unitless
        (this%inhibition_biomass_constant/ &
        (rt_auxvar%immobile(this%biomass_id) + &
        ! A monod inhibition term, unitless
         this%inhibition_biomass_constant)**this%inhibition_biomass_exponent)* &
        (this%inhibition_alcohol_constant/ &
        (this%inhibition_alcohol_constant + &
         rt_auxvar%total(idof_alcohol,iphase)))

  mu_D = this%Cr_food_direct_constant*sum_food*rt_auxvar%total(idof_Cr,iphase)    ! mol/L/s

  Residual(idof_Cr) =      Residual(idof_Cr) + &
                           ! Biological reaction, mol/s
                           this%Cr_decay_constant* &                              ! /s
                           rt_auxvar%immobile(this%biomass_id)* &                 ! mol/m3 bulk
                           rt_auxvar%total(idof_Cr,iphase)/ &                     ! mol/L
                           (this%Cr_inhibition_constant + &                       
                           rt_auxvar%total(idof_Cr,iphase))* &                    ! mol/L
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%Cr_food_Cr_stoic_constant*mu_D* &                 ! mol/L water/s
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol                                  ! L water/m3 bulk

  biomass_residual_delta = &                                                      ! Growth usage, mol/s
                           - mu_B*material_auxvar%volume + &                      ! mol/m3 bulk/s * m3 bulk
                           ! Natural decay, mol/s
                           this%biomass_decay_constant* &                         ! 1/s
                           (rt_auxvar%immobile(this%biomass_id) - &
                            this%biomass_min)* &                                  ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Biocide reaction, mol/s
                           this%biomass_biocide_decay_constant* &                 ! L/mol/s
                           (rt_auxvar%immobile(this%biomass_id) - &
                           this%biomass_min)* &                                   ! mol/m3 bulk
                           rt_auxvar%total(idof_biocide,iphase)* &                ! mol/L
                           material_auxvar%volume                                 ! m3 bulk

  Residual(idof_biomass) = Residual(idof_biomass) + biomass_residual_delta        ! mol/s
  
  mobile_mole_fraction = rt_auxvar%total(idof_food_mobile,iphase)/sum_food
  immobile_mole_fraction = 1 - mobile_mole_fraction


  Residual(idof_food_mobile) = Residual(idof_food_mobile) + &
                           ! Growth usage, mol/s
                           this%food_decay_constant_1* &
                           mobile_mole_fraction* &                                ! dimensionless
                           mu_B*material_auxvar%volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! Direct usage, mol/s
                           this%food_decay_constant_2* &                          ! 1/s
                           mobile_mole_fraction* &                                ! dimensionless
                           rt_auxvar%immobile(this%biomass_id)* &                 ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%Cr_food_food_stoic_constant* &
                           mobile_mole_fraction* &                                ! dimensionless
                           mu_D*material_auxvar%volume* &                         ! mol/L/s * m3 bulk
                           immobile_to_water_vol + &                              ! L water/m3 bulk
                           ! immobilization, mol/s
                           this%food_immobilization_constant* &                   ! 1/s
                           rt_auxvar%total(idof_food_mobile,iphase)* &            ! mol/L
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol - &                              ! L water/m3 bulk
                           ! remobilization, mol/s
                           this%food_remobilization_constant* &                   ! 1/s
                           rt_auxvar%immobile(this%food_immobile_id)* &           ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk


  Residual(idof_food_immobile) = Residual(idof_food_immobile) + &
                           ! Growth usage, mol/s
                           this%food_decay_constant_1* &                          ! unitless
                           immobile_mole_fraction* &                              ! dimensionless
                           mu_B*material_auxvar%volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! Direct usage, mol/s
                           this%food_decay_constant_2* &                          ! 1/s
                           immobile_mole_fraction* &                              ! dimensionless
                           rt_auxvar%immobile(this%biomass_id)* &                 ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%Cr_food_food_stoic_constant* &
                           immobile_mole_fraction* &                              ! dimensionless
                           mu_D*material_auxvar%volume* &                         ! mol/L/s * m3 bulk
                           immobile_to_water_vol - &                              ! L water/m3 bulk
                           ! immobilization, mol/s
                           this%food_immobilization_constant* &                   ! 1/s
                           rt_auxvar%total(idof_food_mobile,iphase)* &            ! mol/L
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol + &                              ! L water/m3 bulk
                           ! remobilization, mol/s
                           this%food_remobilization_constant* &                   ! 1/s
                           rt_auxvar%immobile(this%food_immobile_id)* &           ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk

Residual(idof_biocide) = Residual(idof_biocide) + &
                         this%biocide_decay_constant* &                           ! L/mol/s
                         rt_auxvar%total(idof_biocide,iphase)* &                  ! mol/L
                         rt_auxvar%immobile(this%biomass_id)* &                   ! mol/m3 bulk
                         material_auxvar%volume                                   ! m3 bulk
                         
if (compute_derivative) then
  option%io_buffer = 'Analytical Jacobian calculation for this is a pain.' // &
    ' Use Numerical Jacobian only!!!!'
  call printErrMsg(option)
endif

end subroutine ChromeAlcoholPermReact

! ************************************************************************** !

subroutine ChromeAlcoholPermKineticState(this,rt_auxvar,global_auxvar, &
                                               material_auxvar, &
                                               reaction,option)
  !
  ! Updates the kinetic state for the sandbox
  !
  ! Author: Satish Karra, Scott Hansen, Sachin Pandey, LANL
  ! Date: 09/15/2016
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only: material_auxvar_type

  implicit none

  class(reaction_sandbox_Cr_bio_perm_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: biomass_residual_delta, delta_volfrac
  
  PetscReal :: L_water
  PetscReal :: mu_B, mu_D
  PetscReal :: dmu_dF, dmu_dB, sum_food
  PetscInt :: idof_food_mobile, idof_food_immobile, idof_biomass, idof_Cr
  PetscInt :: idof_alcohol, idof_biomineral, idof_biocide
  PetscReal :: immobile_to_water_vol
  PetscReal :: immobile_mole_fraction, mobile_mole_fraction
  
  idof_food_mobile = this%food_mobile_id
  idof_Cr = this%Cr_id
  idof_alcohol = this%alcohol_id
  idof_biocide = this%biocide_id
  idof_biomass = reaction%offset_immobile + this%biomass_id
  idof_food_immobile = reaction%offset_immobile + this%food_immobile_id
  
  immobile_to_water_vol = &
            material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0                   ! L water/ m3 bulk
  
  sum_food = rt_auxvar%total(idof_food_mobile,iphase) + &
            rt_auxvar%immobile(this%food_immobile_id)/ &
            immobile_to_water_vol                                                ! in mol/L water; Note that food_immobile is divided by porosity*saturation
  
  mu_B = this%biomass_growth_constant*rt_auxvar%immobile(this%biomass_id)* &      ! mol/m3 bulk/s
            ! F monod term, unitless
            (sum_food/(sum_food + this%monod_food_constant))* &
            ! B monod inhibition term, unitless
            (this%inhibition_biomass_constant/ &
            (rt_auxvar%immobile(this%biomass_id) + &
            ! A monod inhibition term, unitless
                this%inhibition_biomass_constant)**this%inhibition_biomass_exponent)* &
            (this%inhibition_alcohol_constant/ &
            (this%inhibition_alcohol_constant + &
                rt_auxvar%total(idof_alcohol,iphase)))
  
  biomass_residual_delta = &                                       ! Growth usage, mol/s
            - mu_B*material_auxvar%volume + &                      ! mol/m3 bulk/s * m3 bulk
            ! Natural decay, mol/s
            this%biomass_decay_constant* &                         ! 1/s
            (rt_auxvar%immobile(this%biomass_id) - &
            this%biomass_min)* &                                   ! mol/m3 bulk
            material_auxvar%volume + &                             ! m3 bulk
            ! Biocide reaction, mol/s
            this%biomass_biocide_decay_constant* &                 ! L/mol/s
            (rt_auxvar%immobile(this%biomass_id) - &
            this%biomass_min)* &                                   ! mol/m3 bulk
            rt_auxvar%total(idof_biocide,iphase)* &                ! mol/L
            material_auxvar%volume                                 ! m3 bulk

  delta_volfrac = &
            - biomass_residual_delta / &                             ! mol/s
            this%biomass_density / &                               ! mol/m3
            material_auxvar%volume * &                             ! m3 bulk
            option%tran_dt                                         ! s

  rt_auxvar%mnrl_volfrac(this%biomineral_id) = rt_auxvar%mnrl_volfrac(this%biomineral_id) + &
                                      delta_volfrac
 
end subroutine ChromeAlcoholPermKineticState

! ************************************************************************** !

subroutine ChromeAlcoholPermDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  implicit none

  class(reaction_sandbox_Cr_bio_perm_type) :: this

! 12. Add code to deallocate contents of the Cr_bio object

end subroutine ChromeAlcoholPermDestroy

end module Reaction_Sandbox_Chrome_Alcohol_Perm_class
