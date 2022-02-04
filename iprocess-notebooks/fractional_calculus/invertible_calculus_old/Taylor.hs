module Taylor (plot, printTaylor, PPlot.range, PPlot.repeatList) where

import qualified PPlot as PPlot

toFunction taylorSeries x = sum [a*x^^n/fromIntegral (PPlot.factorial n) | (a, n) <- zip (take PPlot._taylor_cutoff taylorSeries) [0..] ]

eval taylorSeries xs = [(x, toFunction taylorSeries x) | x<-xs]

plot x functions = PPlot.plot [(eval function x, name) | (function, name) <- functions]
printTaylor taylorSeries = putStrLn $ "Taylor Series: " ++ (init . show . take 10) taylorSeries ++ " ... ]"
