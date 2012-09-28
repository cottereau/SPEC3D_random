module logical_input

! Modified by Gaetano Festa 31/01/2005
! Modified by Elise Delavaud 15/02/2006

type :: Logical_array 

logical :: save_trace,save_snapshots,save_energy,plot_grid,save_restart,run_restart
logical :: run_exec, run_debug, run_echo
logical :: any_source
logical :: super_object, Neumann, super_object_local_present, Neumann_local_present,Save_Surface

end type 

end module logical_input
