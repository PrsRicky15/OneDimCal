program modern_quantum_solver
    use quantum_solver
    use potential
    implicit none
    type(ResonancePoten) :: jol
    jol = create_resonancePoten()
    call real_Hamiltonian(jol)
end program modern_quantum_solver