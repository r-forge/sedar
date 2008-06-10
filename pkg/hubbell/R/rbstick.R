"rbstick" <-
    function(n)
{
    stick <- runif(n-1)
    stick <- c(0,sort(stick),1)
    broke <- stick[-1] - stick[-(n+1)]
    broke <- rev(sort(broke))
    broke
}
