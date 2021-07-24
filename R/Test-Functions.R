#' @export
soma=function(x1,x2)
{
  times4(x1+x2)
}

#' @export
soma2=function(x1,x2)
{
  times3(x1+x2)
}

#' @export
somaC <- function(x, y) .Call(testeC, x, y)
