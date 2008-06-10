"print.landscape" <-
  function(x, ...)
{
  print(table(as.vector(x)))
  invisible(x)
}
