matmult                  package:base                  R Documentation

_M_a_t_r_i_x _M_u_l_t_i_p_l_i_c_a_t_i_o_n

_D_e_s_c_r_i_p_t_i_o_n:

     Multiplies two matrices, if they are conformable.  If one argument
     is a vector, it will be promoted to either a row or column matrix
     to make the two arguments conformable.  If both are vectors it
     will return the inner product (as a matrix).

_U_s_a_g_e:

     x %*% y
     
_A_r_g_u_m_e_n_t_s:

    x, y: numeric or complex matrices or vectors.

_D_e_t_a_i_l_s:

     When a vector is promoted to a matrix, its names are not promoted
     to row or column names, unlike ‘as.matrix’.

     This operator is S4 generic but not S3 generic.  S4 methods need
     to be written for a function of two arguments named ‘x’ and ‘y’.

_V_a_l_u_e:

     A double or complex matrix product.  Use ‘drop’ to remove
     dimensions which have only one level.

_R_e_f_e_r_e_n_c_e_s:

     Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) _The New S
     Language_.  Wadsworth & Brooks/Cole.

_S_e_e _A_l_s_o:

     ‘matrix’, ‘Arithmetic’, ‘diag’.

_E_x_a_m_p_l_e_s:

     x <- 1:4
     (z <- x %*% x)    # scalar ("inner") product (1 x 1 matrix)
     drop(z)             # as scalar
     
     y <- diag(x)
     z <- matrix(1:12, ncol = 3, nrow = 4)
     y %*% z
     y %*% x
     x %*% z
     

