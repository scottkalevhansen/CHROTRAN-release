module Init_Subsurface_Flow_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: InitSubsurfFlowSetupRealization
  
contains

! ************************************************************************** !

subroutine InitSubsurfFlowSetupRealization(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Init_Common_module
  use Material_module
  
  use Flash2_module
  use Mphase_module
  use Immis_module
  use Miscible_module
  use Richards_module
  use TH_module
  use General_module
  use TOilIms_module
  use Condition_Control_module
  use co2_sw_module, only : init_span_wagner
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  
  ! initialize FLOW
  ! set up auxillary variable arrays
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(TH_MODE)
        call THSetup(realization)
      case(RICHARDS_MODE)
        call MaterialSetup(realization%patch%aux%Material%material_parameter, &
                           patch%material_property_array, &
                           patch%characteristic_curves_array, &
                           realization%option)
        call RichardsSetup(realization)
      case(MPH_MODE)
        call init_span_wagner(option)      
        call MphaseSetup(realization)
      case(IMS_MODE)
        call init_span_wagner(option)      
        call ImmisSetup(realization)
      case(MIS_MODE)
        call MiscibleSetup(realization)
      case(FLASH2_MODE)
        call init_span_wagner(option)      
        call Flash2Setup(realization)
      case(G_MODE)
        call MaterialSetup(realization%patch%aux%Material%material_parameter, &
                           patch%material_property_array, &
                           patch%characteristic_curves_array, &
                           realization%option)
        call GeneralSetup(realization)
      case(TOIL_IMS_MODE)
        call MaterialSetup(realization%patch%aux%Material%material_parameter, &
                           patch%material_property_array, &
                           patch%characteristic_curves_array, &
                           realization%option)
        call TOilImsSetup(realization)
    end select
  
    ! assign initial conditionsRealizAssignFlowInitCond
    call CondControlAssignFlowInitCond(realization)

    ! override initial conditions if they are to be read from a file
    if (len_trim(option%initialize_flow_filename) > 1) then
      call InitSubsurfFlowReadInitCond(realization, &
                                       option%initialize_flow_filename)
    endif
  
    select case(option%iflowmode)
      case(TH_MODE)
        call THUpdateAuxVars(realization)
      case(RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case(MPH_MODE)
        call MphaseUpdateAuxVars(realization)
      case(IMS_MODE)
        call ImmisUpdateAuxVars(realization)
      case(MIS_MODE)
        call MiscibleUpdateAuxVars(realization)
      case(FLASH2_MODE)
        call Flash2UpdateAuxVars(realization)
      case(G_MODE)
        !geh: cannot update state during initialization as the guess will be
        !     assigned as the initial conditin if the state changes. therefore,
        !     pass in PETSC_FALSE
        call GeneralUpdateAuxVars(realization,PETSC_FALSE)
      case(TOIL_IMS_MODE)
        call TOilImsUpdateAuxVars(realization)
    end select
  else ! no flow mode specified
    if (len_trim(realization%nonuniform_velocity_filename) > 0) then
#if defined(PETSC_HAVE_HDF5)
      call InitCommonReadVelocityField(realization)
#else
      write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                               &"read HDF5 formatted fluxes in for transport with no flow.")')
      call printErrMsg(option)
#endif
    endif
  endif  
#ifdef WELL_CLASS
  call AllWellsSetup(realization)
#endif
  
end subroutine InitSubsurfFlowSetupRealization

! ************************************************************************** !
#ifdef WELL_CLASS
subroutine AllWellsSetup(realization)
  ! 
  ! Point well auxvars to the domain auxvars te wells belong to
  ! does nothing if well are not defined
  ! 
  ! Author: Paolo Orsini
  ! Date: 06/03/16
  ! 
  use Realization_Subsurface_class
  use Coupler_module
  use Well_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(coupler_type), pointer :: source_sink
  PetscInt :: cpl_idx_start

  source_sink => realization%patch%source_sink_list%first

  cpl_idx_start = 1
  do
    if (.not.associated(source_sink)) exit
    if( associated(source_sink%well) ) then
      !exclude empty wells - not included in well comms
      if(source_sink%connection_set%num_connections > 0) then
        source_sink%well%name = source_sink%name    
        call WellAuxVarSetUp(source_sink%well,source_sink%connection_set, &
                         source_sink%flow_condition,realization%patch%aux, &
                         cpl_idx_start,realization%patch%ss_flow_vol_fluxes, &
                         realization%option)
      end if
    end if
    cpl_idx_start = cpl_idx_start + source_sink%connection_set%num_connections
    source_sink => source_sink%next
  end do

end subroutine AllWellsSetup
#endif

! ************************************************************************** !

subroutine InitSubsurfFlowReadInitCond(realization,filename)
  ! 
  ! Assigns flow initial condition from HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/05/10, 12/04/14
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Discretization_module
  use HDF5_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: local_id, idx, offset
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  if (option%iflowmode /= RICHARDS_MODE) then
    option%io_buffer = 'Reading of flow initial conditions from HDF5 ' // &
                       'file (' // trim(filename) // &
                       'not currently not supported for mode: ' // &

                       trim(option%flowmode)
  endif      

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

      ! assign initial conditions values to domain
    call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)

    ! Pressure for all modes 
    offset = 1
    group_name = ''
    dataset_name = 'Pressure'
    call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                      filename,group_name, &
                                      dataset_name,option%id>0)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id=1, grid%nlmax
      if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
      if (dabs(vec_p(local_id)) < 1.d-40) then
        print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
              ': Potential error - zero pressure in Initial Condition read from file.'
      endif
      idx = (local_id-1)*option%nflowdof + offset
      xx_p(idx) = vec_p(local_id)
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)

    call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
        
    cur_patch => cur_patch%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)  
  call VecCopy(field%flow_xx, field%flow_yy, ierr);CHKERRQ(ierr)

end subroutine InitSubsurfFlowReadInitCond

end module Init_Subsurface_Flow_module
