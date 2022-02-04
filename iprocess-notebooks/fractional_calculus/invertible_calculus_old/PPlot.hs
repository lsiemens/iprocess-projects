module PPlot where

import qualified Graphics.Rendering.Chart.Easy as Easy
import qualified Graphics.Rendering.Chart.Backend.Cairo as Cairo

_taylor_cutoff = 50
_resolution = 100

range :: (Double, Double) -> [Double]
range (start, stop) =
    let dx = (stop - start)/fromIntegral (_resolution - 1)
    in [fromIntegral i*dx + start| i <- [0 .. (_resolution - 1)]]

repeatList = concat . repeat

factorial 0 = 1
factorial n = n*factorial (n - 1)

-- plotting
binding :: (Monad m) => [m a] -> m a
binding [a] = a
binding (a:as) = a >> binding as

plot functions = Easy.toRenderable $ binding [Easy.plot (Easy.line name [points]) | (points, name) <- functions]
