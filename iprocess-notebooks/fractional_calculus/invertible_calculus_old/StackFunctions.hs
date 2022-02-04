module StackFunctions (Differentiable, StackFunction(..), plot, PPlot.range, PPlot.repeatList) where

import qualified PPlot as PPlot

class Differentiable a where
    derivative :: a -> a
    antiderivative :: a -> a

type Function = (Double -> Double)

class Analytic a where
    toFunction :: a -> Function

eval :: (Analytic a) => a -> [Double] -> [(Double, Double)]
eval function xs = [(x, toFunction function x) | x<-xs]

data StackFunction = Constant Double | StackFunction [Double] [Double]

instance Show StackFunction where
    show (StackFunction (a_0:terms) intterms) = "[... " ++ (init.tail.show.reverse.take 4) intterms ++ ",〈" ++ show a_0 ++ "〉," ++ (init.tail.show.take 4) terms ++ " ...]"
    show (Constant a) = "Constant: " ++ show a

instance Num StackFunction where
    (+) (Constant a) (Constant b) = Constant (a + b)
    (+) (Constant a) (StackFunction (b:terms) intterms) = StackFunction (a+b:terms) intterms
    (+) function (Constant a) = Constant a + function 
    (+) (StackFunction termsA inttermsA) (StackFunction termsB inttermsB) = StackFunction [a+b | (a, b) <- zip (termsA ++ repeat 0) (termsB ++ repeat 0)] [a+b | (a, b) <- zip (inttermsA ++ repeat 0) (inttermsB ++ repeat 0)]

    negate (Constant a) = Constant $ negate a
    negate (StackFunction terms intterms) = StackFunction (map negate terms) (map negate intterms)

    abs (Constant a) = Constant $ abs a
    abs (StackFunction terms intterms) = StackFunction (map abs terms) (map abs intterms)

    fromInteger a = Constant $ fromInteger a

    (*) (Constant a) (Constant b) = Constant (a*b)
    (*) (Constant a) (StackFunction terms intterms) = StackFunction (map (a*) terms) (map (a*) intterms)
    (*) function (Constant a) = Constant a * function 
    (*) _ _ = error "Undefined"
    signum _ = error "Undefined"

instance Fractional StackFunction where
    fromRational a = Constant $ fromRational a

    (/) (Constant a) (Constant b) = Constant (a/b)
    (/) (StackFunction terms intterms) (Constant a) = StackFunction (map (*recip a) terms) (map (*recip a) intterms)
    (/) _ _ = error "Undefined"

instance Analytic StackFunction where
    toFunction (StackFunction terms intterms) =
        let function :: StackFunction -> Function
            function (StackFunction terms _) x = sum [a*x^^n/fromIntegral (PPlot.factorial n) | (a, n) <- zip terms [0..] ]
        in function $ StackFunction (take PPlot._taylor_cutoff terms) []

plot x functions = PPlot.plot [(eval function x, name) | (function, name) <- functions]
