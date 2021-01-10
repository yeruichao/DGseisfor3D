module test_mod

    use datatype_mod, only : rkind,matrices,tetmesh_geometry,&
                             vector_array,tensor_array,&
                             reset_vector,reset_tensor,&
                             vector_axpy,tensor_axpy,DiagMM
contains

subroutine gradU_E_diff(pNp,mesh,matrix,U,E,resE)
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    integer :: pNp
    type(vector_array) :: U
    type(tensor_array) :: E,resE
end subroutine gradU_E_diff

end module test_mod
