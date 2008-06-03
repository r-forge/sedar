`plot.chooseCN` <-
    function(x, ...)
{
    NextMethod(x, "plot", coord = attr(x, "xy"), ...)
}
