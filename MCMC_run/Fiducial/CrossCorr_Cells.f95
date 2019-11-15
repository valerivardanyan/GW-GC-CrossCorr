SUBROUTINE C_ells(C_ell_lst, ell_lst, k_unique_lst, eta_unique_lst, bess_int_keta_lst,&
  confH, sigma, deltaCDM, W_kernel, K_kernel, bias, number_ells, k_length, eta_length)
  !This MODULE calculates the C_ells for grav wave backgrounds


  IMPLICIT NONE

  INTEGER :: i, number_ells, ell
  INTEGER, INTENT(IN) :: eta_length, k_length

  DOUBLE PRECISION, INTENT(INOUT) :: C_ell_lst(1:number_ells)
  INTEGER, INTENT(IN) :: ell_lst(1:number_ells)
  DOUBLE PRECISION, INTENT(IN) :: confH(1:k_length, 1:eta_length)
  DOUBLE PRECISION, INTENT(IN) :: sigma(1:k_length, 1:eta_length), deltaCDM(1:k_length, 1:eta_length)
  DOUBLE PRECISION, INTENT(IN) :: W_kernel(1:eta_length), K_kernel(1:eta_length), bias(1:eta_length)

  DOUBLE PRECISION, INTENT(IN) :: k_unique_lst( 1:k_length), eta_unique_lst(1:eta_length)
  DOUBLE PRECISION, INTENT(IN) :: bess_int_keta_lst(1:number_ells, 1:k_length, 1:eta_length)


  DOUBLE PRECISION :: eta_integral_1, eta_integral_2, k_integral
  DOUBLE PRECISION :: k_integral_lst(k_length)
  DOUBLE PRECISION :: W_kernel_tmp, K_kernel_tmp, bias_tmp
  DOUBLE PRECISION :: deltaCDM_tmp, confH_tmp, sigma_tmp, besselInt_tmp

  INTEGER :: k_index, eta_index

  DOUBLE PRECISION :: delta_eta, delta_logk, k, k_next



  DO i = 1, number_ells, 1

    !WRITE(*,*) "Performing computation for l = ", ell_lst(i)

    k_integral = 0.


    DO k_index = 1, k_length, 1

      k = k_unique_lst(k_index)
      k_next = k_unique_lst(k_index + 1)

      eta_integral_1 = 0.
      eta_integral_2 = 0.

      DO eta_index = 2, eta_length, 1

          delta_eta = eta_unique_lst(eta_index) - eta_unique_lst(eta_index - 1)


          W_kernel_tmp = (W_kernel(eta_index) + W_kernel(eta_index - 1))/2.
          K_kernel_tmp = (K_kernel(eta_index) + K_kernel(eta_index - 1))/2.

          bias_tmp = (bias(eta_index) + bias(eta_index - 1))/2.


          deltaCDM_tmp = (deltaCDM(k_index, eta_index) + deltaCDM(k_index, eta_index - 1))/2.
          confH_tmp = (confH(k_index, eta_index) + confH(k_index, eta_index - 1))/2.
          sigma_tmp = (sigma(k_index, eta_index) + sigma(k_index, eta_index - 1))/2.

          besselInt_tmp = -(bess_int_keta_lst(i, k_index, eta_index) - bess_int_keta_lst(i, k_index, eta_index - 1))/k
          eta_integral_1 = eta_integral_1 + W_kernel_tmp*besselInt_tmp*(deltaCDM_tmp - 3.*confH_tmp*sigma_tmp/k/bias_tmp)
          eta_integral_2 = eta_integral_2 + K_kernel_tmp*besselInt_tmp*(deltaCDM_tmp - 3.*confH_tmp*sigma_tmp/k/bias_tmp)


      END DO

      k_integral_lst(k_index) = k_integral + eta_integral_1*eta_integral_2

      !write(*,*) k_integral_lst(k_index)

    END DO

    k_integral = 0.
    DO k_index = 1, k_length - 1, 1

      k = k_unique_lst(k_index)
      k_next = k_unique_lst(k_index + 1)
      delta_logk = DLOG(k_next) - DLOG(k)

      k_integral = k_integral + (k_integral_lst(k_index + 1) + k_integral_lst(k_index))*delta_logk/2.

    END DO


    ell = ell_lst(i)
    C_ell_lst(i) = k_integral


  END DO

END SUBROUTINE C_ells
