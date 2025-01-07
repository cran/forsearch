#' @export
logist3 <-
function (x, a, b, c) 
{
#                                logist3

# VALUE 3-parameter logistic function
#
xx <- exp(-(x - b)/c)
y <- a/(1 + xx)

return(y)
}
