program FUSS_docgen
    use FUSS_Read_Sim_Param, only: Register_Sim_Param
    use FUSS_Read_IO,        only: Register_IO_Fields, Register_Probes
    use FUSS_Read_Numerics,  only: Register_Numerics
    use FUSS_Backend_INI,    only: Load_Ini, Scan_Ini
    use FUSS_Input_Registry
    implicit none

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(1, 'probe-placeholder')  ! Register one probe with placeholder name
    call Register_Numerics(2)

    call reg%generate_markdown('docs/user/registry.md')

end program FUSS_docgen
    
