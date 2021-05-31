#' Teste funcao soma
#' 
#' Essa funcao soma dois numeros
#' 
#' Essa funcao eh apenas para teste do desenvolvimento do pacote R - RobustBetaReg
#' 
#' @param x1 eh um numero
#' @param x2 eh um numero
#' 
#' @return retorna a soma
#' 
#' @export
soma = function(x1,x2)
{
  y=x1+x2
  return(y)
}

#' Testa 
#' 
#' Testa funcao soma positivo
#' 
#' @param x resultado da funcao.
#' @export
teste.soma = function(x)
{
  return(x>0)
}

#' Dados de exemplo
#'
#' Dados para ser utilizado de exemplo
#'   
#' @format .RData
#' 
#' @source mtcars' 
"mtcars2"