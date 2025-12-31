module Data.NonEmpty (NonEmpty, initNE, updNE, getNE) where

-- | A non-empty container for stuff of type a is an item of type a along with a (possibly empty) list of the same type
--
-- Is this already defined somewhere ?
-- data NonEmpty a = NE
--   !a   -- ^ Current state
--   [a] -- ^ Buffer of previous states
--   !Int -- ^ Length of buffer
--   deriving (Eq, Show)

data NonEmpty a = NE !a [a] !Int deriving (Eq, Show)

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

-- | Reconstruct the most recent list of states from a 'NonEmpty' record
getNE :: Int -> NonEmpty a -> Maybe [a]
getNE n (NE sc scs lbuf)
  | lbuf < n - 1 = Nothing
  | lbuf == n = Just (sc : scs)
  | otherwise = Just (sc : init scs)

