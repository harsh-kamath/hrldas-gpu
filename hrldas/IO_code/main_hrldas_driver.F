program Noah_hrldas_driver
! this is the main program to drive HRLDAS-Noah, HRLDAS-NoahMP, and other Land models.

#ifdef Noah1d
! this is used to drive Noah1d
   use module_noah1d_hrldas_driver, only: land_driver_ini, land_driver_exe
#else
! this is used to drive NoahMP
   use module_noahmp_hrldas_driver, only: land_driver_ini, land_driver_exe, land_driver_finish
#ifdef MPP_LAND
   use mpi, only: MPI_Initialized, MPI_Finalized, MPI_Finalize, MPI_SUCCESS
#endif
#endif

   implicit none
   integer :: ITIME,NTIME
#ifndef Noah1d
#ifdef MPP_LAND
   integer :: mpi_error
   logical :: mpi_started, mpi_finished
#endif
#endif
   
   call land_driver_ini(NTIME)
   
   do ITIME = 0, NTIME
       call land_driver_exe(ITIME)
   end do

#ifndef Noah1d
   call land_driver_finish()
#ifdef MPP_LAND
   ! This standalone program owns MPI shutdown. Complete GPU/output cleanup
   ! before the collective finalization; coupled host drivers are unaffected.
   call MPI_Initialized(mpi_started,mpi_error)
   if (mpi_error /= MPI_SUCCESS) error stop 'HRLDAS: MPI_Initialized failed at shutdown'
   if (mpi_started) then
      call MPI_Finalized(mpi_finished,mpi_error)
      if (mpi_error /= MPI_SUCCESS) error stop 'HRLDAS: MPI_Finalized failed at shutdown'
      if (.not. mpi_finished) then
         call MPI_Finalize(mpi_error)
         if (mpi_error /= MPI_SUCCESS) error stop 'HRLDAS: MPI_Finalize failed at shutdown'
      endif
   endif
#endif
#endif
END
