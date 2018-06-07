module Control.Misc (Counter, getCounter, initCounter, incCounter) where








-- | A Counter is a very specialized integer value which always starts at 0 and can only increase by 1
newtype Counter = C {getCounter :: Int} deriving (Eq)
instance Show Counter where show (C i) = show i

initCounter :: Counter
initCounter = C 0

incCounter :: Counter -> Counter
incCounter (C i) = C (i + 1)
