summary.dhglm <-
function (object,...) {
    ans<-list(0)
    ans$FixCoef<-object$FixCoef
    ans$RandCoef<-object$RandCoef
    ans$likelihood<-object$likelihood
    ans$iter<-object$iter
    ans$convergence<-object$convergence
    class(ans) <- "summary.frialtyHL"
    ans
}

