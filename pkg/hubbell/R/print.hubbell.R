"print.hubbell" <-
function (x, ...) 
{
  print(table(x), ...)
  invisible(x)
}
