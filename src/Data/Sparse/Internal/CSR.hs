module Data.Sparse.Internal.CSR where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VM

data CsrMatrix a =
  CsrMatrix { csrVal :: VU.Vector a,
              csrColInd :: VU.Vector Int,
              csrRowPtr :: VU.Vector Int,
              csrNnz :: {-# UNPACK #-} !Int,
              csrNrows :: {-# UNPACK #-} !Int,
              csrNcols :: {-# UNPACK #-} !Int } deriving Eq

