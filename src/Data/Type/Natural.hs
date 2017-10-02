-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Type.Natural
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  BSD3 (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- Type-level naturals inspired by  https://hackage.haskell.org/package/type-natural but available in similar form in many other Haskell packages
--
-----------------------------------------------------------------------------
module Data.Type.Natural where



-- data Zero = Zero deriving (Eq, Show) 
-- data Succ a = Succ a deriving (Eq, Show)
-- type D1 = Succ Zero
-- type D2 = Succ D1
-- type D3 = Succ D2 

-- class Nat a where
--   fromNat :: a -> Int
--   toNat :: Int -> a
-- instance Nat Zero where
--   fromNat _ = 0
--   toNat _ = Zero
-- instance Nat a => Nat (Succ a) where
--   fromNat _ = 1 + fromNat (undefined :: a)
--   toNat n = Succ (toNat (n - 1))

data Nat a = Zero | Succ (Nat a) deriving (Eq, Show)
fromNat :: Num t => Nat t1 -> t
fromNat Zero = 0
fromNat (Succ n) = 1 + fromNat n

toNat 0 = Zero
toNat n = Succ (toNat (n - 1))
-- instance Num (Nat a) where
--   m + n = fromNat m + fromNat n
