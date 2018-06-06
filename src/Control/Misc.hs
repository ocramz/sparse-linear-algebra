module Control.Misc (NonEmpty, initNE, updNE, getNE, Counter, getCounter, initCounter, incCounter) where


-- | A non-empty container for stuff of type a is an item of type a along with a (possibly empty) list of the same type
--
-- Is this already defined somewhere ?
data NonEmpty a = NE {neCurrent :: !a, nePrev :: ![a], nePrevLen :: !Int } deriving (Eq, Show)

initNE :: a -> NonEmpty a
initNE a = NE a [] 0

-- | Update a 'NonEmpty' record.
--
-- 
updNE ::
     Int  -- ^ Window length (NB : /must/ be positive)
  -> a    -- ^ New state
  -> NonEmpty a -- ^ Current
  -> NonEmpty a -- ^ Updated
updNE n s (NE sc scs lbuf) 
  | lbuf < n = NE s (sc : scs) (lbuf + 1)
  | otherwise = NE s (sc : init scs) lbuf

getNE :: Int -> NonEmpty a -> Maybe [a]
getNE n (NE sc scs lbuf)
  | lbuf < n - 1 = Nothing
  | lbuf == n = Just (sc : scs)
  | otherwise = Just (sc : init scs)

buf3 :: a -> NonEmpty a -> NonEmpty a
buf3 = updNE 3    





-- | A Counter is a very specialized integer value which always starts at 0 and can only increase by 1
newtype Counter = C {getCounter :: Int} deriving (Eq)
instance Show Counter where show (C i) = show i

initCounter :: Counter
initCounter = C 0

incCounter :: Counter -> Counter
incCounter (C i) = C (i + 1)
