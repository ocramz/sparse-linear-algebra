module Numeric.LinearAlgebra.Core.Sparse where

import Data.List (unfoldr)
import Control.Monad.State.Strict

import Numeric.LinearAlgebra.Core.Class





rotMtx m n ixrow ixcol angle = undefined
  where
    -- mat = smFromList (m, n)


rotMtx0 m n irow jcol angle = undefined where
  arrx = [0 .. m-1]
  arry = [0.. n-1]
  c = cos angle
  s = sin angle
  -- go coords
  --   | coords == (irow, irow) = c
  --   | coords == (irow, jcol) = s
  --   | coords == (jcol, jcol) = c
  --   | coords == (jcol, irow) = - s
  --   | otherwise = coords
  


data Entry a = E { cx :: Int, cy :: Int, elem :: a}


step m angle ii jj (i, j)
  | i < m - 1 || j < m - 1 =
        if i /= ii && j /= jj
        then Just (E (i+1) (j+1) 1, (i + 1, j + 1))
        else undefined -- if i == ii then 
  | otherwise = Nothing

