program modern_quantum_solver
    use potential
    use quantum_solver
    implicit none
    type(ResonancePoten) :: jol
    jol = create_resonancePoten()
    call real_Hamiltonian(jol)
end program modern_quantum_solver