crPlotsHGLM <-
function (model, terms = ~., layout = NULL, ask, main, ...) 
{
    terms <- if (is.character(terms)) 
        paste("~", terms)
    else terms
    vform <- update(formula(model), terms)
    if (any(is.na(match(all.vars(vform), all.vars(formula(model)))))) 
        stop("Only predictors in the formula can be plotted.")
    mf <- attr(model.frame(model), "terms")
    terms <- attr(mf, "term.labels")
    vterms <- attr(terms(vform), "term.labels")
    if (any(attr(terms(model), "order") > 1)) {
        stop("C+R plots not available for models with interactions.")
    }
    nt <- length(vterms)
    if (nt == 0) 
        stop("No plots specified")
    if (missing(main)) 
        main <- if (nt == 1) 
            "cubic spline smoothing"
        else "cubic spline smoothing"
    if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
        if (is.null(layout)) {
            layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 
                2), c(2, 2), c(3, 2), c(3, 2), c(3, 3), c(3, 
                3), c(3, 3))
        }
        ask <- if (missing(ask) || is.null(ask)) 
            prod(layout) < nt
        else ask
        op <- par(mfrow = layout, ask = ask, no.readonly = TRUE, 
            oma = c(0, 0, 1.5, 0), mar = c(5, 4, 1, 2) + 0.1)
        on.exit(par(op))
    }
    if (!is.null(class(model$na.action)) && class(model$na.action) == 
        "exclude") 
        class(model$na.action) <- "omit"
    for (term in vterms) res<-crPlot(model, term, ...)
    # mtext(side = 3, outer = TRUE, main, cex = 1.2)
    # invisible(0)
    return(res)
}
