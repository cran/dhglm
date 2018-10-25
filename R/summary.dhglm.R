summary.dhglm <-
function (object) {
    print("Estimates from the model(mu)")
    print(object$beta_coeff,4)
    print("Estimates for logarithm of lambda=var(u_mu)")
    print(object$lambda_coeff,4)
    print("Estimates from the model(phi)")
    print(object$phi_coeff)
    print("========== Likelihood Function Values and Condition AIC ==========")
    print(object$likeli_coeff)
}
